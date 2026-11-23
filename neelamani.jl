### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ af3c2bae-e340-11f1-bae1-bd2a249053d4
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


# ╔═╡ 526f4c8a-21ff-41d8-831d-9dcdea7aab2a
begin
	qs = [50, 60, 70, 75,80, 85, 90, 95]
	jpegs = Dict(q => Float64.(rawview(channelview(load("./data/grey_exp/train01_q$(q).png")))) for q in qs)
	png = Float64.(rawview(channelview(load("./data/grey_exp/train01_grey.png"))))
end;

# ╔═╡ be8c4c62-51b9-4527-b567-e01bd7ff3936
# Chunk up image into 8x8 blocks
blocks(x) = [x[i:i+7, j:j+7] for i in 1:8:8(size(x,1)÷8), j in 1:8:8(size(x,2)÷8)];

# ╔═╡ e1f008dc-f1ea-4b91-a27f-d1826aaf98de
show(image) = Gray.(image ./ maximum(image))

# ╔═╡ 1f09e3c2-be35-4072-977a-e0fcbb4b4303
show(blocks(jpegs[50])[100])

# ╔═╡ 3556d243-0c72-444e-92ab-398c9a1ff919
begin
	Y(image) = [dct(B) for B in blocks(image)]
	c01(image) = [y[1, 2] for y in Y(image)]
end

# ╔═╡ 6fc2a78b-689b-4499-8ae1-79cc4efc45ab
Y(jpegs[50])

# ╔═╡ 087d6998-a5fb-42d5-b2f4-2c3d56976b1e
show(c01(jpegs[50]))

# ╔═╡ a3f784c2-8829-4c89-a72d-a2ab37b6ed8e
begin
	  c_hist(c, vals) = [count(==(v), round.(Int, c)) for v in vals] ./ length(c)
	vals = -60:60
	p = plot(vals, c_hist(c01(png), vals), label="uncompressed", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p, vals, c_hist(c01(jpegs[q]), vals), label="q=$(q)")
	end
	p
end

# ╔═╡ e83664e1-8344-4de5-9428-5288525b5e67

begin
	p2 = plot(vals, max.(c_hist(c01(png), vals), 1e-6), label="uncompressed", title="Histogram of Y(0,1)", yscale=:log10)
	for q in qs 
		plot!(p2, vals, max.(c_hist(c01(jpegs[q]), vals), 1e-4), label="q=$(q)", yscale=:log10)
	end
	p2
end

# ╔═╡ 58dfa365-7932-488e-b32a-908ac03e8240
begin
	p3 = plot(abs.(rfft(c_hist(c01(png), vals))), label="uncompressed", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p3, abs.(rfft(c_hist(c01(jpegs[q]), vals))), label="q=$(q)")
	end
	p3
end

# ╔═╡ ef2f0225-eb11-47c4-9479-1e2c101e5236
begin
	p4 = plot([0], label="uncompressed", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p4, abs.(rfft(c_hist(c01(jpegs[q]), vals))) - abs.(rfft(c_hist(c01(png), vals))), label="q=$(q)")
	end
	p4
end

# ╔═╡ Cell order:
# ╠═af3c2bae-e340-11f1-bae1-bd2a249053d4
# ╠═526f4c8a-21ff-41d8-831d-9dcdea7aab2a
# ╠═be8c4c62-51b9-4527-b567-e01bd7ff3936
# ╠═e1f008dc-f1ea-4b91-a27f-d1826aaf98de
# ╠═1f09e3c2-be35-4072-977a-e0fcbb4b4303
# ╠═3556d243-0c72-444e-92ab-398c9a1ff919
# ╠═6fc2a78b-689b-4499-8ae1-79cc4efc45ab
# ╠═087d6998-a5fb-42d5-b2f4-2c3d56976b1e
# ╠═a3f784c2-8829-4c89-a72d-a2ab37b6ed8e
# ╠═e83664e1-8344-4de5-9428-5288525b5e67
# ╠═58dfa365-7932-488e-b32a-908ac03e8240
# ╠═ef2f0225-eb11-47c4-9479-1e2c101e5236
