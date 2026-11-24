

#import "@preview/wordometer:0.1.5": word-count, total-words 

#show: word-count
Word count: #total-words / 4000
= DSP Project

#outline()
== Introduction
== Greyscale Images
=== Estimating Compression from Cross-Boundary Differences

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
  caption: [Fan–de Queiroz cross-boundary difference for each greyscale test image.],
) <tab:cross-boundary>


=== Estimating Quantisation Factor $Q$

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

=== Single-Compression vs Double-Compression

#figure(
  text(size: 7pt, table(
    columns: 10,
    align: (left,) + (right,) * 8 + (left,),
    table.header(
      table.cell(rowspan: 2, align: horizon)[Image],
      table.cell(colspan: 8, align: center)[DCT coefficient $(u, v)$],
      table.cell(rowspan: 2, align: horizon)[Compression verdict],
      [$(0, 1)$], [$(0, 2)$], [$(1, 0)$], [$(1, 1)$], [$(1, 2)$], [$(2, 0)$], [$(2, 1)$], [$(2, 2)$],
    ),
    [grey01], [--], [--], [--], [0.0002 (31)], [0.0066 (13)], [0.0211 (13)], [--], [--], [double-compressed],
    [grey02], [0.0020 (30)], [0.0025 (27)], [0.0033 (17)], [0.0032 (23)], [0.0010 (24)], [0.0013 (27)], [0.0049 (31)], [0.0031 (26)], [uncompressed],
    [grey03], [0.0014 (27)], [0.0000 (29)], [0.0002 (29)], [--], [--], [--], [--], [--], [single-compressed],
    [grey04], [--], [--], [--], [0.0000 (31)], [--], [--], [--], [--], [single-compressed],
    [grey05], [0.0043 (25)], [0.0016 (26)], [0.0027 (17)], [0.0037 (28)], [0.0009 (18)], [0.0029 (24)], [0.0027 (30)], [0.0007 (29)], [uncompressed],
    [grey06], [0.0032 (14)], [0.0027 (31)], [0.0020 (25)], [0.0009 (31)], [--], [--], [0.0002 (27)], [--], [uncompressed],
    [grey07], [0.0069 (15)], [0.0009 (31)], [0.0029 (31)], [0.0016 (28)], [--], [0.0012 (28)], [0.0031 (31)], [--], [uncompressed],
    [grey08], [--], [--], [--], [--], [--], [--], [--], [--], [single-compressed],
    [grey09], [0.0048 (10)], [0.0035 (27)], [0.0039 (21)], [0.0049 (20)], [0.0015 (20)], [0.0028 (29)], [0.0025 (31)], [--], [single-compressed],
    [grey10], [0.1506 (19)], [0.0644 (22)], [0.1574 (19)], [0.1004 (19)], [0.0480 (31)], [0.0839 (31)], [0.0433 (25)], [--], [double-compressed],
    [grey11], [--], [0.0086 (31)], [--], [--], [--], [--], [--], [0.0112 (30)], [double-compressed],
    [grey12], [0.0375 (31)], [--], [0.0024 (31)], [--], [0.0034 (23)], [0.0317 (24)], [0.0011 (28)], [--], [double-compressed],
    [grey13], [--], [0.0001 (31)], [--], [--], [--], [--], [--], [--], [single-compressed],
    [grey14], [--], [0.0005 (31)], [0.0007 (29)], [--], [--], [0.0014 (31)], [--], [--], [uncompressed],
    [grey15], [0.0032 (25)], [0.0041 (30)], [0.0046 (31)], [0.0024 (18)], [0.0031 (17)], [0.0025 (22)], [0.0029 (29)], [0.0319 (31)], [single-compressed],
  )),
  caption: [Double-compression peak score for each greyscale test image (rows) and DCT coefficient $(u, v)$ (columns), with the index of the peak in the histogram spectrum in parentheses; -- means no peak above the spectral envelope. Images ruled uncompressed in @tab:cross-boundary are labelled uncompressed; the remaining images are labelled double-compressed if their requantised histograms show periodic notches (@fig:double-hist), and single-compressed otherwise.],
) <tab:double-q>

#figure(
  image("double_compression_histogram.svg", width: 100%),
    caption: [Histograms of the requantised coefficients $c_(u v) / q_2$ for each greyscale test image. Periodic notches indicate double compression.],
) <fig:double-hist>


