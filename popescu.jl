### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ 32e2bf36-e749-11f1-bcb7-a3a853c4441d
begin
    import Pkg
    Pkg.activate()  
    using DSP, FFTW, Plots, Statistics
    using DSP
	using Plots
	using FFTW
	using Statistics
	using LinearAlgebra
	using ImageShow 
	using Images
end


# ╔═╡ eefe1cb6-90e3-4253-8640-017fe282c27c
begin
	qs = [50, 60, 70, 75,80, 85, 90, 95]
	ns = 1:5
	load_grey(path) = Float64.(rawview(channelview(load(path))))
	name(n) = "train$(lpad(n, 2, '0'))"
	# Images are loaded lazily from disk on each access. Preloading all
	# 5*(8+64+1) = 365 images as Float64 (~24 MB each) needs ~9 GB of RAM.
	struct LazyImages{F}
		path::F
	end
	Base.getindex(L::LazyImages, k...) = load_grey(L.path(k...))
	single_jpegs = LazyImages((n, q) -> "./data/grey_exp/$(name(n))_q$(q).png")
	double_jpegs = LazyImages((n, q1, q2) -> "./data/grey_double/$(name(n))_q$(q1)_q$(q2).png")
	pngs = LazyImages(n -> "./data/grey_exp/$(name(n))_grey.png")
end;

# ╔═╡ ea971183-d0c9-4148-9b31-2a09dd50e954
# Chunk up image into 8x8 blocks
blocks(x) = [x[i:i+7, j:j+7] for i in 1:8:8(size(x,1)÷8), j in 1:8:8(size(x,2)÷8)];

# ╔═╡ 539b8fd3-f82d-442f-8fb9-757757e31b3a
show(image) = Gray.(image ./ maximum(image))

# ╔═╡ b5eb1354-44c9-45a3-b9ad-3625a3498b6e
begin
	Y(image) = [dct(B) for B in blocks(image)]
	# Orthonormal 1-D DCT-II basis vector (same normalisation as FFTW's dct / JPEG)
	dct_basis(u) = [(u == 0 ? sqrt(1/8) : sqrt(2/8)) * cos((2x + 1) * u * π / 16) for x in 0:7]
	# Single coefficient (u,v) of every 8x8 block, without computing the full DCT
	function cuv(image, u, v)
		bu, bv = dct_basis(u), dct_basis(v)
		[dot(bu, view(image, i:i+7, j:j+7), bv)
		 for i in 1:8:8(size(image,1)÷8), j in 1:8:8(size(image,2)÷8)]
	end
	c01(image) = cuv(image, 0, 1)
end


# ╔═╡ 7192fe57-1a4f-44cf-ae26-ecd123c20719

function pf_spectrum(image, u, v, q2; K = 30)
    k = round.(Int, cuv(image, u, v) ./ q2)
    h = [count(==(n), k) for n in -K:K] ./ length(k)
    # abs.(rfft(h .- mean(h)))
    abs.(rfft(h))
end

# ╔═╡ af4bdf3f-f1cd-442c-82d2-80dd732db512
begin
	# Discrete Laplacian (two-sided geometric) fitted by maximum likelihood to the
	# quantised values k, evaluated on -K:K. Same normalisation as h.
	function laplace_fit(k, K)
		p0 = mean(k .== 0)
		r  = 1 - 1 / mean(abs.(k[k .!= 0]))
		[n == 0 ? p0 : (1 - p0) * (1 - r) * r^(abs(n) - 1) / 2 for n in -K:K]
	end

	# Spectrum of the histogram after removing the fitted Laplacian envelope
	function residual_spectrum(image, u, v, q2; K = 30)
		k = round.(Int, cuv(image, u, v) ./ q2)
		h = [count(==(n), k) for n in -K:K] ./ length(k)
		abs.(rfft(h .- laplace_fit(k, K)))
	end
end


# ╔═╡ 15dfa605-2ae7-4c68-861d-c27810806929
begin
    const T_LUM = [16 11 10 16 24 40 51 61;
                   12 12 14 19 26 58 60 55;
                   14 13 16 24 40 57 69 56;
                   14 17 22 29 51 87 80 62;
                   18 22 37 56 68 109 103 77;
                   24 35 55 64 81 104 113 92;
                   49 64 78 87 103 121 120 101;
                   72 92 95 98 112 100 103 99]
    
    function qtable(Q)
        s = Q < 50 ? 5000 ÷ Q : 200 - 2Q
        clamp.((T_LUM .* s .+ 50) .÷ 100, 1, 255)
    end
    qstep(Q, u, v) = qtable(Q)[u+1, v+1]
end

# ╔═╡ 5e6157fa-c8f3-41e2-b760-cf1e93aae98a
begin
		  c_hist(c, vals) = [count(==(v), round.(Int, c)) for v in vals] ./ length(c)
		vals = -60:60
end

# ╔═╡ fdac42f3-0c2a-445b-b03c-7b07e3fb3730

begin
	q_level_1 = 75
	n_img = 1
	p2 = plot(vals, max.(c_hist(c01(single_jpegs[n_img, q_level_1]), vals), 1e-4), label="q=$(q_level_1)", title="Histogram of Y(0,1)", yscale=:log10)
	for q in qs 
		plot!(p2, vals, max.(c_hist(c01(double_jpegs[n_img, q_level_1, q]), vals), 1e-4), label="q=$((q_level_1, q))", yscale=:log10)
	end
	p2
end

# ╔═╡ a95e67ab-a811-4c56-a2c8-2d72edf736e5

begin
	p3 = plot(abs.(rfft(c_hist(c01(single_jpegs[n_img, q_level_1]), vals))), label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p3, abs.(rfft(c_hist(c01(double_jpegs[n_img, q_level_1, q]), vals))), label="q=$((q_level_1, q))")
	end
	p3
end

# ╔═╡ f2d6c622-5e08-49f2-b5af-b3224ac57b21
#pf_spectrum(image, u, v, q2; K = 30)
begin
	p6 = plot(pf_spectrum(single_jpegs[n_img, q_level_1], 0, 1, qstep(q_level_1, 0, 1)), label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p6, pf_spectrum(double_jpegs[n_img, q_level_1, q], 0, 1, qstep(q, 0, 1)), label="q=$((q_level_1, q))")
	end
	p6
end

# ╔═╡ 714ab274-125e-4662-8182-364094eaed29
# Overlay of all single- and double-compressed images.
# Colour: R = 0 single / 1 double, B = first quality, G = second quality (0 if single)
begin
	qnorm(q) = (q - minimum(qs)) / (maximum(qs) - minimum(qs))
	p7 = plot(title="PF spectrum of Y(0,1), all images", xlabel="frequency bin", legend=false, xlims=(0, Inf))
	for n in ns, q1 in qs, q2 in qs
		plot!(p7, pf_spectrum(double_jpegs[n, q1, q2], 0, 1, qstep(q2, 0, 1)),
			color=RGB(1, qnorm(q2), qnorm(q1)), alpha=0.4)
	end
	for n in ns, q in qs
		plot!(p7, pf_spectrum(single_jpegs[n, q], 0, 1, qstep(q, 0, 1)),
			color=RGB(0, 0, qnorm(q)), alpha=0.6, linewidth=2)
	end
	p7
end

# ╔═╡ 31730b1b-d0b0-4ca0-ad2e-fe0e1ced66e0
# Same overlay and colours as p7, but with the fitted Laplacian removed first
begin
	p8 = plot(title="Residual spectrum of Y(0,1), all images", xlabel="frequency bin", legend=false, xlims=(0, Inf))
	for n in ns, q1 in qs, q2 in qs
		plot!(p8, residual_spectrum(double_jpegs[n, q1, q2], 0, 1, qstep(q2, 0, 1)),
			color=RGB(1, qnorm(q2), qnorm(q1)), alpha=0.4)
	end
	for n in ns, q in qs
		plot!(p8, residual_spectrum(single_jpegs[n, q], 0, 1, qstep(q, 0, 1)),
			color=RGB(0, 0, qnorm(q)), alpha=0.6, linewidth=2)
	end
	p8
end

# ╔═╡ b6f47eef-d354-4cfb-811b-69ab91682894

begin
	p4 = plot([0], label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p4, abs.(rfft(c_hist(c01(double_jpegs[n_img, q_level_1, q]), vals))) .- abs.(rfft(c_hist(c01(single_jpegs[n_img, q]), vals))), label="q=$((q_level_1, q))")
	end
	p4
end

# ╔═╡ 1a707304-9f9d-41ed-b6e0-a0e48636d67e
begin
	p5 = plot([0], label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p5, abs.(rfft(c_hist(c01(double_jpegs[n_img, q, q_level_1]), vals))) .- abs.(rfft(c_hist(c01(single_jpegs[n_img, q_level_1]), vals))), label="q=$((q, q_level_1))")
	end
	p5
end

# ╔═╡ e43ad705-836c-4e0c-a8cf-045c80061e49
# Histograms for several AC coefficients, one panel per
# original image, with every single/double-compressed derivative overlaid.
# Colour: R = 0 single / 1 double, B = (last) JPEG quality, G = AC coefficient
let K = 15
	ac_coefs = [(0, 1), (1, 0), (1, 1), (0, 2), (2, 0)]
	cnorm(i) = (i - 1) / max(length(ac_coefs) - 1, 1)
	function add_hist!(p, img, q, col)
		for (ci, c) in enumerate(ac_coefs)
			k = round.(Int, cuv(img, c...) ./ qstep(q, c...))
			h = [count(==(j), k) for j in -K:K] ./ length(k)
			plot!(p, -K:K, max.(h, 1e-5); color = col(cnorm(ci)))
		end
	end
	newplot(t) = plot(title = t, yscale = :log10, legend = false)
	# Row n: [single-compressed | double-compressed], both derived from trainNN
	ps = mapreduce(vcat, ns) do n
		ps_single = newplot("$(name(n)) single")
		for q in qs
			add_hist!(ps_single, single_jpegs[n, q], q, g -> RGB(0, g, qnorm(q)))
		end
		ps_double = newplot("$(name(n)) double")
		for q1 in qs, q2 in qs
			add_hist!(ps_double, double_jpegs[n, q1, q2], q2, g -> RGB(1, g, qnorm(q2)))
		end
		[ps_single, ps_double]
	end
	plot(ps...; layout = (length(ns), 2), size = (1200, 300 * length(ns)),
	     link = :y,
	     plot_title = "Histograms of round(c_uv / q2), AC coefficients")
end

# ╔═╡ Cell order:
# ╠═32e2bf36-e749-11f1-bcb7-a3a853c4441d
# ╠═eefe1cb6-90e3-4253-8640-017fe282c27c
# ╠═ea971183-d0c9-4148-9b31-2a09dd50e954
# ╠═539b8fd3-f82d-442f-8fb9-757757e31b3a
# ╠═b5eb1354-44c9-45a3-b9ad-3625a3498b6e
# ╠═7192fe57-1a4f-44cf-ae26-ecd123c20719
# ╠═af4bdf3f-f1cd-442c-82d2-80dd732db512
# ╠═15dfa605-2ae7-4c68-861d-c27810806929
# ╠═5e6157fa-c8f3-41e2-b760-cf1e93aae98a
# ╠═fdac42f3-0c2a-445b-b03c-7b07e3fb3730
# ╠═a95e67ab-a811-4c56-a2c8-2d72edf736e5
# ╠═f2d6c622-5e08-49f2-b5af-b3224ac57b21
# ╠═714ab274-125e-4662-8182-364094eaed29
# ╠═31730b1b-d0b0-4ca0-ad2e-fe0e1ced66e0
# ╠═b6f47eef-d354-4cfb-811b-69ab91682894
# ╠═1a707304-9f9d-41ed-b6e0-a0e48636d67e
# ╠═e43ad705-836c-4e0c-a8cf-045c80061e49
