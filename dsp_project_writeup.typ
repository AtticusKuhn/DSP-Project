#import "@preview/codelst:2.0.2": sourcecode

#import "@preview/wordometer:0.1.5": word-count, total-words 
#show: word-count
#let argmax = $op("arg max", limits: #true)$
Word count: #total-words / 4000

CRSid: ak2518

Author: Atticus Kuhn

= DSP Project: Identifying JPEG Compression in Images

#outline()
== Introduction
I read and implemented the techniques in Julia from the works of
Fan @Fan, Neelamani @Neelamani, and Popescu @Popescu
on forensic analysis of JPEG-compressed images.

Each of the existing works presents a different technique
for JPEG analysis:
- Fan presents a technique for detecting JPEG compression from cross-block differences, and an MLE for the quantisation factor $Q$.
- Neelamani explains the histogram patterns of JPEG-compressed images.
- Popescu explains a pattern that characterises double-compressed images.
- Lewis @Lewis provides an algorithm for exact JPEG recompression. However, I did not implement this algorithm.


I calibrated the techniques on a synthetic dataset of images with known compression histories, and then
I applied those techniques on a testing dataset of real-world images of unknown compression (15 coloured and 15
greyscale). I classified each image into one of three categories:
1. Uncompressed
2. Single-compressed
3. Double-compressed

== Greyscale Images
=== Creating a Synthetic Dataset

In order to calibrate the existing techniques, and to demonstrate the reliability of the techniques, I
created a synthetic dataset of images with known compression histories.

`./data/grey_exp/train$NN_q$Q.png` is training image $N$, compressed to quality level $Q$.

`./data/grey_double/train$NN_q$Q1_q$Q2.png` is training image $N$, first compressed to quality level $Q_1$, and then compressed to quality level $Q_2$.

The dataset was generated from the five training images provided with the assignment using 
@code:bash.

#figure(
    sourcecode[```bash
in_dir=$1
out_dir=$2
qualities=(50 60 70 75 80 85 90 95)
mkdir -p "$out_dir"
for f in "$in_dir"/train*.png; do
  name=$(basename "$f" .png)
  grey="$out_dir/${name}_grey.png"
  pngtopnm "$f" | ppmtopgm | pnmtopng > "$grey"
  for q in "${qualities[@]}"; do
    pngtopnm "$grey" | cjpeg -quality "$q" | djpeg | pnmtopng > "$out_dir/${name}_q$(printf %02d "$q").png"
  done
  echo "done: $name"
done```],caption:[Bash script for generating the synthetic dataset]
)<code:bash>


=== Estimating Compression from Cross-Boundary Differences

Fan and de Queiroz demonstrate that in a JPEG-compressed image, pixel differences
across block boundaries are statistically different from pixel differences
within blocks. They define a cross-boundary metric $Z$, which I implemented in @code:fan_Z_K.


#figure(
    sourcecode[```julia
Z(image) = abs.(conv(image, [1.0 -1; -1 1]))
Z1_Z2(image) = let 
	Z_png = Z(image)
	H, W = size(Z_png)
	p = mod1.(1:H, 8)            
	q = mod1.(1:W, 8)            
	(Z_png[p .== 4, q .== 4],Z_png[p .== 8, q .== 8])
end;
zhist(z) = [count(==(v), z) for v in 0:510] ./ length(z)
# In Fan and de Queiroz, \[K = \sum_n | H_I(n) - H_{II}(n)|\]
FQ_K(image) = let (z1, z2) = Z1_Z2(image)
    sum(abs.(zhist(z1) .- zhist(z2)))
end
    ```],caption:[Julia code for computing Fan's $Z$ metric]
)<code:fan_Z_K>

I calibrated the following decision rule:
images with $K > K_("thresh") = 0.03$ are ruled to be
JPEG-compressed.
This threshold is robust on the training set, as evidenced by 
@tab:confusion.
With this decision procedure, I produced the verdicts in 
@tab:cross-boundary.


#figure(
  image("k_vs_jpeg_quality.svg", width: 100%),
    caption: [$K$ vs JPEG quality],
) <fig:k-jpeg>

#figure(
  table(
    columns: 3,
    align: (left, center, center),
    table.header(
      [], [Actually JPEG], [Actually uncompressed],
    ),
    [Ruled JPEG], [40], [1],
    [Ruled uncompressed], [0], [4],
  ),
  caption: [Confusion matrix on synthetic training dataset.],
) <tab:confusion>




#figure(
  table(
    columns: 3,
    align: (center, right, left),
    table.header[Image][Cross-boundary difference][Compression ruling],
    [grey01], [0.6740], [JPEG compressed],
    [grey02], [0.0171], [Uncompressed],
    [grey03], [0.7080], [JPEG compressed],
    [grey04], [0.9267], [JPEG compressed],
    [grey05], [0.0205], [Uncompressed],
    [grey06], [0.0248], [Uncompressed],
    [grey07], [0.0175], [Uncompressed],
    [grey08], [1.0445], [JPEG compressed],
    [grey09], [0.5891], [JPEG compressed],
    [grey10], [0.9444], [JPEG compressed],
    [grey11], [1.0502], [JPEG compressed],
    [grey12], [0.8734], [JPEG compressed],
    [grey13], [1.0047], [JPEG compressed],
    [grey14], [0.0111], [Uncompressed],
    [grey15], [0.4974], [JPEG compressed],
  ),
    caption: [Fan & de Queiroz cross-boundary difference for each greyscale test image.],
) <tab:cross-boundary>


=== Estimating Quantisation Factor $Q$

Fan and de Queiroz prove a formula for the MLE of $hat(q)(u,v)$, which is the quantisation
factor in JPEG compression.

Their method uses the following formula, which I implemented in @code:qhat.

$ r(i,s,m,n,q) = floor.l (Y(i)[m, n][s]) / q floor.r $
$ G(x) = e^(-6 x^2) $
$ L(c,m, n, q) = sum_(y in c) log(sum_(k = -K)^K G(y - floor.l y / q floor.r q + k q)) + |c| log(q)  $
$ hat(q)(i, m, n) = argmax_(q) L("coeffs"(i, m, n), m, n, q) $

#figure(
    sourcecode[```julia
begin
	r(image,s ,m, n, q) = round(Y(image)[m,n][s] / q)
	D(m) = m == 0 || m == 4 ? 2 : m == 2 || m == 6 ? 2cos(pi/4) : 2cos(pi/8)cos(pi/8)
	B(m, n) = D(m)D(n)
	G(x) = exp( - 6 * x^2)
	function L(c, m, n, q)                     
        b = B(m, n)
        K = div(b, q)
        sum(log(sum(G(y - round(y/q)*q - k*q) for k in -K:K))  
            for y in c) + length(c) * log(q)                                    # + N log q
	end
	qhat(image, m, n) = (c = coeffs(image, m, n); argmax(q -> L(c, m, n, q), 1:60))
end
    ```], caption: [Julia implementation of estimating $hat(q)$.]) <code:qhat>

As we can observe from @fig:Q-vs-Qhat, the estimated $hat(Q)$ is a good approximation for the true $Q$.
#figure(
  image("Q_vs_Qhat.svg", width: 100%),
    caption: [$Q$ vs $hat(Q)$],
) <fig:Q-vs-Qhat>


In @tab:q-estimate, I estimate $Q$ (or $Q_2$ in the case of a double-compressed image) for each testing image.
#figure(
  text(size: 8pt, table(
    columns: 12,
    align: (left,) + (right,) * 11,
    table.header(
      table.cell(rowspan: 2, align: horizon)[Image],
      table.cell(rowspan: 2, align: horizon)[$hat(Q)$],
      table.cell(colspan: 9, align: center)[Per-coefficient $Q$ estimate],
      table.cell(rowspan: 2, align: horizon)[Spread],
      [1], [2], [3], [4], [5], [6], [7], [8], [9],
    ),
    [grey01], [84.3], [84.4], [81.8], [83.3], [83.3], [85.0], [85.7], [85.7], [84.6], [84.4], [3.90],
    [grey02\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [grey03], [58.5], [59.4], [59.1], [58.3], [58.3], [60.0], [57.1], [57.1], [57.7], [59.4], [2.86],
    [grey04], [65.3], [65.6], [63.6], [66.7], [66.7], [65.0], [64.3], [64.3], [65.4], [65.6], [3.03],
    [grey05\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [grey06\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [grey07\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [grey08], [80.1], [81.2], [81.8], [79.2], [79.2], [80.0], [78.6], [78.6], [80.8], [81.2], [3.25],
    [grey09], [84.3], [84.4], [81.8], [83.3], [83.3], [85.0], [85.7], [85.7], [84.6], [84.4], [3.90],
    [grey10], [78.8], [78.1], [77.3], [79.2], [79.2], [80.0], [78.6], [78.6], [80.8], [78.1], [3.50],
    [grey11], [87.7], [90.6], [81.8], [83.3], [83.3], [90.0], [89.3], [89.3], [88.5], [90.6], [8.81],
    [grey12], [75.4], [75.0], [77.3], [75.0], [75.0], [75.0], [75.0], [75.0], [76.9], [75.0], [2.27],
    [grey13], [53.8], [53.1], [54.5], [54.2], [54.2], [55.0], [53.6], [53.6], [53.8], [53.1], [1.88],
    [grey14\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [grey15], [84.7], [84.4], [86.4], [83.3], [83.3], [85.0], [85.7], [85.7], [84.6], [84.4], [3.03],
  )),
    caption: [Estimated JPEG quality factor $Q$ for each greyscale test image. Images marked \* were ruled uncompressed.],
) <tab:q-estimate>

=== Single vs Double Compression

Popescu gives a characterisation of certain features associated with double-compressed JPEGs, as opposed
to single-compressed JPEGs.

Popescu claims that the Fourier transform of the histogram of the coefficients of a JPEG will smoothly decay in the form
of a Laplacian. However, if there is double compression, the double-quantisation artefacts will appear as
spikes in the Laplacian envelope. These spikes are visible in @fig:d-vs-s-hist.

// #figure(
//     image("pf_spectrum.png", width: 100%),
//     caption: [PF spectrum of single vs double compressed images],
// ) <fig:Q-vs-Qhat>

#figure(
    image("double_vs_single_histogram.png", width: 100%),
    caption: [Histograms of single-compressed vs double-compressed images],
) <fig:d-vs-s-hist>


I analysed the histograms of the test images by visual inspection. At first, I attempted to develop
a "spikiness/smoothness" metric to quantitatively discriminate between single and
double compression; however, visual inspection proved more reliable.



// #figure(
//   text(size: 7pt, table(
//     columns: 10,
//     align: (left,) + (right,) * 8 + (left,),
//     table.header(
//       table.cell(rowspan: 2, align: horizon)[Image],
//       table.cell(colspan: 8, align: center)[DCT coefficient $(u, v)$],
//       table.cell(rowspan: 2, align: horizon)[Compression verdict],
//       [$(0, 1)$], [$(0, 2)$], [$(1, 0)$], [$(1, 1)$], [$(1, 2)$], [$(2, 0)$], [$(2, 1)$], [$(2, 2)$],
//     ),
//     [grey01], [--], [--], [--], [0.0002 (31)], [0.0066 (13)], [0.0211 (13)], [--], [--], [double-compressed],
//     [grey02], [0.0020 (30)], [0.0025 (27)], [0.0033 (17)], [0.0032 (23)], [0.0010 (24)], [0.0013 (27)], [0.0049 (31)], [0.0031 (26)], [uncompressed],
//     [grey03], [0.0014 (27)], [0.0000 (29)], [0.0002 (29)], [--], [--], [--], [--], [--], [single-compressed],
//     [grey04], [--], [--], [--], [0.0000 (31)], [--], [--], [--], [--], [single-compressed],
//     [grey05], [0.0043 (25)], [0.0016 (26)], [0.0027 (17)], [0.0037 (28)], [0.0009 (18)], [0.0029 (24)], [0.0027 (30)], [0.0007 (29)], [uncompressed],
//     [grey06], [0.0032 (14)], [0.0027 (31)], [0.0020 (25)], [0.0009 (31)], [--], [--], [0.0002 (27)], [--], [uncompressed],
//     [grey07], [0.0069 (15)], [0.0009 (31)], [0.0029 (31)], [0.0016 (28)], [--], [0.0012 (28)], [0.0031 (31)], [--], [uncompressed],
//     [grey08], [--], [--], [--], [--], [--], [--], [--], [--], [single-compressed],
//     [grey09], [0.0048 (10)], [0.0035 (27)], [0.0039 (21)], [0.0049 (20)], [0.0015 (20)], [0.0028 (29)], [0.0025 (31)], [--], [single-compressed],
//     [grey10], [0.1506 (19)], [0.0644 (22)], [0.1574 (19)], [0.1004 (19)], [0.0480 (31)], [0.0839 (31)], [0.0433 (25)], [--], [double-compressed],
//     [grey11], [--], [0.0086 (31)], [--], [--], [--], [--], [--], [0.0112 (30)], [double-compressed],
//     [grey12], [0.0375 (31)], [--], [0.0024 (31)], [--], [0.0034 (23)], [0.0317 (24)], [0.0011 (28)], [--], [double-compressed],
//     [grey13], [--], [0.0001 (31)], [--], [--], [--], [--], [--], [--], [single-compressed],
//     [grey14], [--], [0.0005 (31)], [0.0007 (29)], [--], [--], [0.0014 (31)], [--], [--], [uncompressed],
//     [grey15], [0.0032 (25)], [0.0041 (30)], [0.0046 (31)], [0.0024 (18)], [0.0031 (17)], [0.0025 (22)], [0.0029 (29)], [0.0319 (31)], [single-compressed],
//   )),
//   caption: [Double-compression peak score for each greyscale test image (rows) and DCT coefficient $(u, v)$ (columns), with the index of the peak in the histogram spectrum in parentheses.],
// ) <tab:double-q>

#figure(
  image("double_compression_histogram.svg", width: 100%),
    caption: [Histograms of the requantised coefficients $c_(u v) / q_2$ for each greyscale test image. Periodic notches indicate double compression.],
) <fig:double-hist>

=== `rjpeg` Recompression

`rjpeg`#footnote[#link("https://micrological.appspot.com/cl-web/rjpeg-0.8.tar.gz")] is a tool developed by Lewis for JPEG recompression. I did not write my own implementation of it, but rather
I downloaded and compiled the existing implementation in order to sanity-check my other results. Except in cases
where `rjpeg` times out, it agrees with my other compressed/uncompressed classifications, and with my estimated $Q$
quality factor to within $plus.minus 1$.

#figure(
  text(size: 8pt, table(
    columns: 5,
    align: (left, right, right, right, right),
    table.header(
      table.cell(rowspan: 2, align: horizon)[Image],
      table.cell(rowspan: 2, align: horizon)[`rjpeg` QF],
      table.cell(colspan: 3, align: center)[Y blocks (of 46,128)],
      [Exact], [Ambiguous], [Infeasible],
    ),
    [grey01], [84], [40,525], [5,535], [68],
    [grey02], [--], [0], [0], [46,128],
    [grey03], [58], [45,609], [0], [519],
    [grey04], [65], [46,105], [0], [23],
    [grey05], [--], [0], [0], [46,128],
    [grey06], [--], [0], [0], [46,128],
    [grey07], [--], [0], [0], [46,128],
    [grey08], [80], [43,641], [1,869], [618],
    [grey09], [84], [45,959], [0], [169],
    [grey10], [79], [45,970], [1], [157],
    [grey11], table.cell(colspan: 4, align: center)[timed out after 15 min],
    [grey12], [76], [46,080], [5], [43],
    [grey13], [54], [46,084], [0], [44],
    [grey14], [--], [0], [0], [46,128],
    [grey15], [85], [45,453], [0], [675],
  )),
    caption: [`rjpeg` results for each greyscale test image.],
) <tab:rjpeg>

Looking at the original JPEG gives us a more accurate
understanding of the DCT coefficients, as can be seen in
@fig:rjpeg_pf.

#figure(
  image("rjpeg_pf_spectrum.svg", width: 100%),
    caption: [FFT of coefficients of recompressed images],
) <fig:rjpeg_pf>


== Colour Images

I converted the `./dataset/colourNN.png` images into greyscale (via
RGB $arrow.r$ YCbCr coordinate conversion) and
ran the same greyscale analysis on them. Thus, I will
not repeat the same explanation as in the greyscale
section.

#figure(
  table(
    columns: 3,
    align: (center, right, left),
    table.header[Image][Cross-boundary difference][Compression ruling],
    [colour01], [0.8676], [JPEG compressed],
    [colour02], [0.9847], [JPEG compressed],
    [colour03], [0.0191], [Uncompressed],
    [colour04], [0.5543], [JPEG compressed],
    [colour05], [0.0182], [Uncompressed],
    [colour06], [1.0526], [JPEG compressed],
    [colour07], [0.0207], [Uncompressed],
    [colour08], [1.1144], [JPEG compressed],
    [colour09], [0.6459], [JPEG compressed],
    [colour10], [0.8355], [JPEG compressed],
    [colour11], [0.0210], [Uncompressed],
    [colour12], [0.6469], [JPEG compressed],
    [colour13], [1.1507], [JPEG compressed],
    [colour14], [0.0244], [Uncompressed],
    [colour15], [1.2634], [JPEG compressed],
  ),
    caption: [Fan & de Queiroz cross-boundary difference for each colour test image (converted to greyscale).],
) <tab:cross-boundary-colour>

#figure(
  text(size: 8pt, table(
    columns: 12,
    align: (left,) + (right,) * 11,
    table.header(
      table.cell(rowspan: 2, align: horizon)[Image],
      table.cell(rowspan: 2, align: horizon)[$hat(Q)$],
      table.cell(colspan: 9, align: center)[Per-coefficient $Q$ estimate],
      table.cell(rowspan: 2, align: horizon)[Spread],
      [1], [2], [3], [4], [5], [6], [7], [8], [9],
    ),
    [colour01], [50.0], [50.0], [50.0], [50.0], [50.0], [50.0], [50.0], [50.0], [50.0], [50.0], [0.00],
    [colour02], [88.1], [87.5], [86.4], [87.5], [87.5], [90.0], [89.3], [89.3], [88.5], [87.5], [3.64],
    [colour03\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [colour04], [84.3], [84.4], [81.8], [83.3], [83.3], [85.0], [85.7], [85.7], [84.6], [84.4], [3.90],
    [colour05\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [colour06], [63.6], [90.6], [77.3], [54.2], [54.2], [55.0], [53.6], [53.6], [76.9], [53.1], [37.50],
    [colour07\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [colour08], [69.9], [68.8], [68.2], [70.8], [70.8], [70.0], [71.4], [71.4], [69.2], [68.8], [3.25],
    [colour09], [63.6], [62.5], [63.6], [62.5], [62.5], [65.0], [64.3], [64.3], [65.4], [62.5], [2.88],
    [colour10], [75.4], [75.0], [77.3], [75.0], [75.0], [75.0], [75.0], [75.0], [76.9], [75.0], [2.27],
    [colour11\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [colour12], [92.8], [93.8], [90.9], [91.7], [91.7], [95.0], [92.9], [92.9], [92.3], [93.8], [4.09],
    [colour13], [64.8], [96.9], [59.1], [58.3], [58.3], [60.0], [60.7], [60.7], [61.5], [59.4], [38.54],
    [colour14\*], [96.2], [96.9], [95.5], [95.8], [95.8], [95.0], [96.4], [96.4], [96.2], [96.9], [1.88],
    [colour15], [69.9], [68.8], [68.2], [70.8], [70.8], [70.0], [71.4], [71.4], [69.2], [68.8], [3.25],
  )),
    caption: [Estimated JPEG quality factor $Q$ for each colour test image (converted to greyscale). Images marked \* were ruled uncompressed.],
) <tab:q-estimate-colour>


#figure(
    image("double_compression_histogram_colour.svg", width: 100%),
    caption: [Histograms of the requantised coefficients for each colour test image (converted to greyscale).],
) <fig:hist-coloured>


== Final Conclusions
#figure(
  table(
    columns: 2,
    align: (center, left),
    table.header[Image][Compression ruling],
    [grey01], [Double compressed],
    [grey02],  [Uncompressed],
    [grey03], [Single compressed],
    [grey04],  [Single compressed],
    [grey05],  [Uncompressed],
    [grey06],  [Uncompressed],
    [grey07],  [Uncompressed],
    [grey08],  [Single compressed],
    [grey09],  [Single compressed],
    [grey10], [Double compressed],
    [grey11],  [Double compressed],
    [grey12],  [Double compressed],
    [grey13],  [Single compressed],
    [grey14],  [Uncompressed],
    [grey15],  [Double compressed],
    [colour01], [Single compressed],
    [colour02],  [Double compressed],
    [colour03], [Uncompressed],
    [colour04],  [Single compressed],
    [colour05],  [Uncompressed],
    [colour06],  [Double compressed],
    [colour07],  [Uncompressed],
    [colour08],  [Double compressed],
    [colour09],  [Single compressed],
    [colour10], [Single compressed],
    [colour11],  [Uncompressed],
    [colour12],  [Single compressed],
    [colour13],  [Double compressed],
    [colour14],  [Uncompressed],
    [colour15],  [Single compressed],
  ),
    caption: [Final compression-history verdict for each test image.],
) <tab:final>


#bibliography("bib.bib")
