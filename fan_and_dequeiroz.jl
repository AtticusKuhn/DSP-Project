### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ 1ee3c82e-e335-11f1-10fc-d95d9f4d3297
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



# ╔═╡ 9f83cb50-ce7e-4bb4-9dfc-4a9a34f7f5e5
begin
	qs = [50, 60, 70, 75,80, 85, 90, 95]
	jpegs = Dict(q => Float64.(rawview(channelview(load("./data/grey_exp/train01_q$(q).png")))) for q in qs)
	png = Float64.(rawview(channelview(load("./data/grey_exp/train01_grey.png"))))
end;

# ╔═╡ 7b279050-1d54-4a78-a312-b006f351f24c
show(image) = Gray.(image ./ maximum(image))

# ╔═╡ 63a09d23-4a22-4546-8064-a4332163c713
show(png)

# ╔═╡ b1c89fe0-6afb-432b-9d5b-8d36450edd31
# In the Fan and De Queiroz Paper, Z is defined as Z = |A - B - C + D|
Z(image) = abs.(conv(image, [1.0 -1; -1 1]));

# ╔═╡ af07c7cb-e647-4747-95dd-acdec110bd2a
show(Z(png))


# ╔═╡ 18229cc0-2cbf-4719-8300-805e33ed086c
begin
	Z_png = Z(png)
	H, W = size(Z_png)
	p = mod1.(1:H, 8)            # position of each row of Z within its block
	q = mod1.(1:W, 8)            # position of each column of Z within its block
	
	Z1 = Z_png[p .== 4, q .== 4]     # within-block centre  (same as Z[4:8:end, 4:8:end])
	Z2 = Z_png[p .== 8, q .== 8]     # across the grid corner (same as Z[8:8:end, 8:8:end])
	
	# To see where these samples are, colour them on top of the image:
	
	# img = RGB.(Gray.(png ./ 255))
	# img[4:8:end, 4:8:end] .= RGB(0, 1, 0)    # Z′ positions: green
	# img[8:8:end, 8:8:end] .= RGB(1, 0, 0)    # Z″ positions: red
	# img
	# img[1:64, 1:64]                          # zoom in on an 8×8 grid of blocks
end;

# ╔═╡ d3d7bf64-83dd-430b-83a4-45ca58038b05
begin
	img = RGB.(Gray.(png ./ 255))
	img[4:8:end, 4:8:end] .= RGB(0, 1, 0)    # Z′ positions: green
	img[8:8:end, 8:8:end] .= RGB(1, 0, 0)    # Z″ positions: red
	img
end

# ╔═╡ 5fc26df7-4916-4c4b-b939-03e30b5e9a5d
Z1_Z2(image) = let 
	Z_png = Z(image)
	H, W = size(Z_png)
	p = mod1.(1:H, 8)            # position of each row of Z within its block
	q = mod1.(1:W, 8)            # position of each column of Z within its block
	
	Z1 = Z_png[p .== 4, q .== 4]     # within-block centre  (same as Z[4:8:end, 4:8:end])
	Z2 = Z_png[p .== 8, q .== 8]     # across the grid corner (same as Z[8:8:end, 8:8:end])
	(Z1, Z2)
end;

# ╔═╡ c3674974-fc7c-4e2f-9fa8-766315c4c7e3
begin
        # normalised histogram over the integer values 0..510
        zhist(z) = [count(==(v), z) for v in 0:510] ./ length(z)

        # In Fan and De Quieroz, \[K = \sum_n | H_I(n) - H_{II}(n)|\]
        FQ_K(image) = let (z1, z2) = Z1_Z2(image)
                sum(abs.(zhist(z1) .- zhist(z2)))
        end

        function plot_FQ(image; xmax=30)
                z1, z2 = Z1_Z2(image)
                plot(0:510, zhist(z1), seriestype=:steppost, label="H_I (within block)")
                plot!(0:510, zhist(z2), seriestype=:steppost, label="H_II (across grid)")
                plot!(xlims=(0, xmax), xlabel="Z", ylabel="relative frequency",
                      title="K = $(round(FQ_K(image), digits=3))")
        end
end

# ╔═╡ ade7d33c-8b81-4033-9f4b-bf561c015d5b
plot_FQ(png)

# ╔═╡ ee6e487d-f135-4347-a7df-c223b248d46a
plot_FQ(jpegs[85])

# ╔═╡ e9f07e01-e31c-46b7-a65b-9efe7ca917db
let
        p = plot(xlabel="Z", ylabel="|H_I − H_II|", xlims=(0, 5))
        d(img) = (z = Z1_Z2(img); (abs.(zhist(z[1]) .- zhist(z[2]))))
        plot!(p, 0:510, d(png), lw=2, color=:black, label="uncompressed")
        for q in qs
                plot!(p, 0:510, d(jpegs[q]), label="q=$q")
        end
        p
end

# ╔═╡ f6aeb3c3-32b1-4cc6-936d-63457f7f7a88
let
        p = plot(xlabel="Z", ylabel="cumulative |H_I − H_II|", xlims=(0, 5))
        d(img) = (z = Z1_Z2(img); cumsum(abs.(zhist(z[1]) .- zhist(z[2]))))
        plot!(p, 0:510, d(png), lw=2, color=:black, label="uncompressed")
        for q in qs
                plot!(p, 0:510, d(jpegs[q]), label="q=$q")
        end
        p
end

# ╔═╡ 7d9808af-2d95-4d25-8022-31aa864f35c1
# Chunk up image into 8x8 blocks
blocks(x) = [x[i:i+7, j:j+7] for i in 1:8:8(size(x,1)÷8), j in 1:8:8(size(x,2)÷8)];

# ╔═╡ f9e3caf2-54ce-40fd-92ea-4849ec9601ab


# ╔═╡ Cell order:
# ╠═1ee3c82e-e335-11f1-10fc-d95d9f4d3297
# ╠═9f83cb50-ce7e-4bb4-9dfc-4a9a34f7f5e5
# ╠═63a09d23-4a22-4546-8064-a4332163c713
# ╠═7b279050-1d54-4a78-a312-b006f351f24c
# ╠═b1c89fe0-6afb-432b-9d5b-8d36450edd31
# ╠═af07c7cb-e647-4747-95dd-acdec110bd2a
# ╠═18229cc0-2cbf-4719-8300-805e33ed086c
# ╠═d3d7bf64-83dd-430b-83a4-45ca58038b05
# ╠═5fc26df7-4916-4c4b-b939-03e30b5e9a5d
# ╠═c3674974-fc7c-4e2f-9fa8-766315c4c7e3
# ╠═ade7d33c-8b81-4033-9f4b-bf561c015d5b
# ╠═ee6e487d-f135-4347-a7df-c223b248d46a
# ╠═e9f07e01-e31c-46b7-a65b-9efe7ca917db
# ╠═f6aeb3c3-32b1-4cc6-936d-63457f7f7a88
# ╠═7d9808af-2d95-4d25-8022-31aa864f35c1
# ╠═f9e3caf2-54ce-40fd-92ea-4849ec9601ab
