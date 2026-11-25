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
	train_ids = 1:5
	loadgrey(path) = Float64.(rawview(channelview(load(path))))
	# jpegs[i, q]: training image i compressed at quality q
	jpegs = Dict((i, q) => loadgrey("./data/grey_exp/train$(lpad(i, 2, '0'))_q$(q).png")
	             for i in train_ids, q in qs)
	# pngs[i]: uncompressed training image i
	pngs = [loadgrey("./data/grey_exp/train$(lpad(i, 2, '0'))_grey.png") for i in train_ids]
end;

# ╔═╡ 7b279050-1d54-4a78-a312-b006f351f24c
show(image) = Gray.(image ./ maximum(image))

# ╔═╡ 63a09d23-4a22-4546-8064-a4332163c713
show(pngs[1])

# ╔═╡ b1c89fe0-6afb-432b-9d5b-8d36450edd31
# In the Fan and De Queiroz Paper, Z is defined as Z = |A - B - C + D|
# Z(image) = abs.(conv(image, [1.0 -1; -1 1]));
Z(x) = abs.(x[1:end-1, 1:end-1] .- x[2:end, 1:end-1] .- x[1:end-1, 2:end] .+ x[2:end, 2:end])

# ╔═╡ af07c7cb-e647-4747-95dd-acdec110bd2a
show(Z(pngs[1]))


# ╔═╡ 62ecbcee-db6e-495a-a8a1-088a67943ad6
# mod1.(1:size(Z(pngs[1])), 8)

# ╔═╡ 18229cc0-2cbf-4719-8300-805e33ed086c
begin
	Z_png = Z(pngs[1])
	H, W = size(Z_png)
	p = mod1.(1:H, 8)          
	q = mod1.(1:W, 8)          
	
	Z1 = Z_png[p .== 4, q .== 4]    
	Z2 = Z_png[p .== 8, q .== 8]     
	
	# To see where these samples are, colour them on top of the image:
	
	# img = RGB.(Gray.(pngs[1] ./ 255))
	# img[4:8:end, 4:8:end] .= RGB(0, 1, 0)    
	# img[8:8:end, 8:8:end] .= RGB(1, 0, 0)   
	# img
	# img[1:64, 1:64]                       
end;

# ╔═╡ d3d7bf64-83dd-430b-83a4-45ca58038b05
begin
	img = RGB.(Gray.(pngs[1] ./ 255))
	img[4:8:end, 4:8:end] .= RGB(0, 1, 0)    # Z′ positions: green
	img[8:8:end, 8:8:end] .= RGB(1, 0, 0)    # Z″ positions: red
	img
end

# ╔═╡ 5fc26df7-4916-4c4b-b939-03e30b5e9a5d
Z1_Z2(image) = let 
	Z_png = Z(image)
	H, W = size(Z_png)
	p = mod1.(1:H, 8)            
	q = mod1.(1:W, 8)            
	
	Z1 = Z_png[p .== 4, q .== 4]      
	Z2 = Z_png[p .== 8, q .== 8]   
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
plot_FQ(pngs[1])

# ╔═╡ ee6e487d-f135-4347-a7df-c223b248d46a
plot_FQ(jpegs[1, 85])

# ╔═╡ e9f07e01-e31c-46b7-a65b-9efe7ca917db
let
        p = plot(xlabel="Z", ylabel="|H_I − H_II|", xlims=(0, 5))
        d(img) = (z = Z1_Z2(img); (abs.(zhist(z[1]) .- zhist(z[2]))))
        plot!(p, 0:510, d(pngs[1]), lw=2, color=:black, label="uncompressed")
        for q in qs
                plot!(p, 0:510, d(jpegs[1, q]), label="q=$q")
        end
        p
end

# ╔═╡ f6aeb3c3-32b1-4cc6-936d-63457f7f7a88
let
        p = plot(xlabel="Z", ylabel="cumulative |H_I − H_II|", xlims=(0, 5))
        d(img) = (z = Z1_Z2(img); cumsum(abs.(zhist(z[1]) .- zhist(z[2]))))
        plot!(p, 0:510, d(pngs[1]), lw=2, color=:black, label="uncompressed")
        for q in qs
                plot!(p, 0:510, d(jpegs[1, q]), label="q=$q")
        end
        p
end

# ╔═╡ 7d9808af-2d95-4d25-8022-31aa864f35c1
# Chunk up image into 8x8 blocks
blocks(x) = [x[i:i+7, j:j+7] for i in 1:8:8(size(x,1)÷8), j in 1:8:8(size(x,2)÷8)];

# ╔═╡ f9e3caf2-54ce-40fd-92ea-4849ec9601ab
md"
# MLE for estimating `q`
"

# ╔═╡ 14485c16-9f3c-4ad0-9a69-9a2f5484e817
md"
Equations taken from Fan and De Quieroz.
```math
Y^*(m, n)= Y(m, n) + noise
```

```math
r_s(m, n, q) = round[Y^*_s(m, n) / q]
```

```math
\hat{q}(m, n) = argmax_q L(q) =  argmax_q { \sum_s \log { \sum_k G(Y^*_s(m, n) − r_s(m, n, q) q − k q) } + N \log q } 
```

```math
L(m, n, q) =   { \sum_s \log { \sum_k G(Y^*_s(m, n) − r_s(m, n, q) q − k q) } + N \log q }
```
```math
B(m,n) = D(m)D(n)
```

```math
D(m) = \text{some nonsense}
```


s = block index. 

k = offset to neighbouring multiples of q. 

```math
G(x) = e^{−6x^2}/Z \qquad|x| ≤ B(m,n)
```
"

# ╔═╡ 1447446c-d959-4b4b-ad17-d48f526b785b
Y(image) = [dct(B) for B in blocks(image)]

# ╔═╡ ff31f989-b442-4488-8619-f2b0614bcf48
r(image,s ,m, n, q) = round(Y(image)[m,n][s] / q)

# ╔═╡ 91af9506-df33-4be8-94d2-dcfc7d136e6b
begin
	# L = log likelhood
	D(m) = m == 0 || m == 2 ? 2 : m == 2 || m == 6 ? 2cos(pi/4) : 2cos(pi/8)cos(pi/8)
	B(m, n) = D(m)D(n)
	G(x) = exp( - 6 * x^2)
	# L(image, m, n, q) = let 
	# hx = div(size(image, 1), 8)
	# hy= div(size(image, 2), 8)
	# N = hx * hy
	# sum(
	# 	log(
	# 		sum(
	# 			G(Y(image)[s][m,n] - r(image, (sx, sy), m ,n,q) * q - k * q)
	# 			+ N * log(q)
	# 			for k in -div(B(m, n),q ):(div(B(m, n),q))
	# 		)
	# 	)
	# 	for sx in 1:hx 
	# 	for sy in 1:hy
	# )
	# end
	function L(c, m, n, q)                     
        b = B(m, n)
        K = div(b, q)
        sum(log(sum(G(y - round(y/q)*q - k*q) for k in -K:K))  # Σ_s Σ_k
            for y in c) + length(c) * log(q)                                    # + N log q
	end
end

# ╔═╡ 488021c6-9576-4d8a-8dd8-f2405c02f254
begin
	  coeffs(image, m, n) = [dct(blk .- 128)[m+1, n+1] for blk in blocks(image)]
	qhat(image, m, n) = (c = coeffs(image, m, n); argmax(q -> L(c, m, n, q), 1:60))
end

# ╔═╡ cd385154-9eea-4005-bbeb-a2fa1f524dec
qhat(jpegs[1, 75], 0, 1)

# ╔═╡ 4ea661d0-095b-46a6-a227-7b0fb58b7237
qhat(jpegs[1, 50], 0, 1)

# ╔═╡ dbae0c1c-364f-43d8-8744-8509fa8b3fd7
let
        plt = plot(xlabel="Q", ylabel="K",
                   title="K vs JPEG quality", legend=:topright)
        for i in train_ids
                xs = [qs; 100]
                ys = [[FQ_K(jpegs[i, q]) for q in qs]; FQ_K(pngs[i])]
                scatter!(plt, xs, ys, color=i, label="train$(lpad(i, 2, '0'))")
        end
        plt
end

# ╔═╡ d5a2eff2-09a8-49d5-b713-b0ea97e97ef4
let
        freqs = [(0,0), (0,1), (1,0), (1,1), (0,2), (2,0), (1,2), (2,1), (2,2)]
        nI, nF = length(train_ids), length(freqs)
        red(i)  = (i - 1) / max(nI - 1, 1)  
        blue(k) = (k - 1) / max(nF - 1, 1)   

        # IJG standard luminance table, and the inverse of libjpeg's quality scaling
        T_LUM = [16  11  10  16  24  40  51  61
                 12  12  14  19  26  58  60  55
                 14  13  16  24  40  57  69  56
                 14  17  22  29  51  87  80  62
                 18  22  37  56  68 109 103  77
                 24  35  55  64  81 104 113  92
                 49  64  78  87 103 121 120 101
                 72  92  95  98 112 100 103  99]
        alpha_to_Q(α) = α <= 1 ? 100 - 50α : 50 / α

        xs, ys, cols = Float64[], Float64[], RGB{Float64}[]
        for i in train_ids
                for (Q, im) in [[(Qj, jpegs[i, Qj]) for Qj in qs]; (100, pngs[i])]
                        Ys = [dct(blk .- 128) for blk in blocks(im)]
                        for (k, (u, v)) in enumerate(freqs)
                                c = [Yb[u+1, v+1] for Yb in Ys]
                                qh = argmax(q -> L(c, u, v, q), 1:60)          # q̂(u,v), the step size
                                push!(xs, Q + 0.4 * (k - (nF + 1) / 2))
                                push!(ys, alpha_to_Q(qh / T_LUM[u+1, v+1]))   # Q̂ from this coefficient
                                push!(cols, RGB(red(i), 0, blue(k)))
                        end
                end
        end

        plt = scatter(xs, ys; markercolor = cols, markerstrokewidth = 0, label = false,
                      xlabel = "true Q (100 = uncompressed)", ylabel = raw"$\hat{Q}$",
                      title = "Estimated vs true JPEG quality",
                      legend = :outerright, size = (900, 550))
        plot!(plt, [40, 100], [40, 100]; color = :grey, ls = :dash, label = "y = x")
        for i in train_ids
                scatter!(plt, [NaN], [NaN]; markercolor = RGB(red(i), 0, 0), markerstrokewidth = 0,
                         label = "train$(lpad(i, 2, '0'))")
        end
        for (k, f) in enumerate(freqs)
                scatter!(plt, [NaN], [NaN]; markercolor = RGB(0, 0, blue(k)), markerstrokewidth = 0,
                         label = "(u,v) = $f")
        end
        plt
end

# ╔═╡ Cell order:
# ╠═1ee3c82e-e335-11f1-10fc-d95d9f4d3297
# ╠═9f83cb50-ce7e-4bb4-9dfc-4a9a34f7f5e5
# ╠═63a09d23-4a22-4546-8064-a4332163c713
# ╠═7b279050-1d54-4a78-a312-b006f351f24c
# ╠═b1c89fe0-6afb-432b-9d5b-8d36450edd31
# ╠═af07c7cb-e647-4747-95dd-acdec110bd2a
# ╠═62ecbcee-db6e-495a-a8a1-088a67943ad6
# ╠═18229cc0-2cbf-4719-8300-805e33ed086c
# ╠═d3d7bf64-83dd-430b-83a4-45ca58038b05
# ╠═5fc26df7-4916-4c4b-b939-03e30b5e9a5d
# ╠═c3674974-fc7c-4e2f-9fa8-766315c4c7e3
# ╠═ade7d33c-8b81-4033-9f4b-bf561c015d5b
# ╠═ee6e487d-f135-4347-a7df-c223b248d46a
# ╠═e9f07e01-e31c-46b7-a65b-9efe7ca917db
# ╠═f6aeb3c3-32b1-4cc6-936d-63457f7f7a88
# ╠═7d9808af-2d95-4d25-8022-31aa864f35c1
# ╟─f9e3caf2-54ce-40fd-92ea-4849ec9601ab
# ╠═14485c16-9f3c-4ad0-9a69-9a2f5484e817
# ╠═1447446c-d959-4b4b-ad17-d48f526b785b
# ╠═ff31f989-b442-4488-8619-f2b0614bcf48
# ╠═91af9506-df33-4be8-94d2-dcfc7d136e6b
# ╠═488021c6-9576-4d8a-8dd8-f2405c02f254
# ╠═cd385154-9eea-4005-bbeb-a2fa1f524dec
# ╠═4ea661d0-095b-46a6-a227-7b0fb58b7237
# ╠═dbae0c1c-364f-43d8-8744-8509fa8b3fd7
# ╠═d5a2eff2-09a8-49d5-b713-b0ea97e97ef4
