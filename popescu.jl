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
	single_jpegs = Dict(q => Float64.(rawview(channelview(load("./data/grey_exp/train01_q$(q).png")))) for q in qs)
	double_jpegs = Dict((q1, q2) => Float64.(rawview(channelview(load("./data/grey_double/train01_q$(q1)_q$(q2).png")))) for q1 in qs for q2 in qs)
	png = Float64.(rawview(channelview(load("./data/grey_exp/train01_grey.png"))))
end;

# ╔═╡ ea971183-d0c9-4148-9b31-2a09dd50e954
# Chunk up image into 8x8 blocks
blocks(x) = [x[i:i+7, j:j+7] for i in 1:8:8(size(x,1)÷8), j in 1:8:8(size(x,2)÷8)];

# ╔═╡ 539b8fd3-f82d-442f-8fb9-757757e31b3a
show(image) = Gray.(image ./ maximum(image))

# ╔═╡ b5eb1354-44c9-45a3-b9ad-3625a3498b6e
begin
	Y(image) = [dct(B) for B in blocks(image)]
	c01(image) = [y[1, 2] for y in Y(image)]
	  cuv(image, u, v) = [y[u+1, v+1] for y in Y(image)]
end


# ╔═╡ 7192fe57-1a4f-44cf-ae26-ecd123c20719

function pf_spectrum(image, u, v, q2; K = 30)
    k = round.(Int, cuv(image, u, v) ./ q2)
    h = [count(==(n), k) for n in -K:K] ./ length(k)
    abs.(rfft(h .- mean(h)))
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
	p2 = plot(vals, max.(c_hist(c01(single_jpegs[q_level_1]), vals), 1e-4), label="q=$(q_level_1)", title="Histogram of Y(0,1)", yscale=:log10)
	for q in qs 
		plot!(p2, vals, max.(c_hist(c01(double_jpegs[(q_level_1, q)]), vals), 1e-4), label="q=$((q_level_1, q))", yscale=:log10)
	end
	p2
end

# ╔═╡ a95e67ab-a811-4c56-a2c8-2d72edf736e5

begin
	p3 = plot(abs.(rfft(c_hist(c01(single_jpegs[q_level_1]), vals))), label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p3, abs.(rfft(c_hist(c01(double_jpegs[(q_level_1, q)]), vals))), label="q=$((q_level_1, q))")
	end
	p3
end

# ╔═╡ f2d6c622-5e08-49f2-b5af-b3224ac57b21
#pf_spectrum(image, u, v, q2; K = 30)
begin
	p6 = plot(pf_spectrum(single_jpegs[q_level_1], 0, 1, qstep(q_level_1, 0, 1)), label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p6, pf_spectrum(double_jpegs[(q_level_1, q)], 0, 1, qstep(q, 0, 1)), label="q=$((q_level_1, q))")
	end
	p6
end

# ╔═╡ b6f47eef-d354-4cfb-811b-69ab91682894

begin
	p4 = plot([0], label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p4, abs.(rfft(c_hist(c01(double_jpegs[(q_level_1, q)]), vals))) .- abs.(rfft(c_hist(c01(single_jpegs[q]), vals))), label="q=$((q_level_1, q))")
	end
	p4
end

# ╔═╡ 1a707304-9f9d-41ed-b6e0-a0e48636d67e
begin
	p5 = plot([0], label="q=$(q_level_1)", title="Histogram of Y(0,1)")
	for q in qs 
		plot!(p5, abs.(rfft(c_hist(c01(double_jpegs[(q, q_level_1)]), vals))) .- abs.(rfft(c_hist(c01(single_jpegs[q_level_1]), vals))), label="q=$((q, q_level_1))")
	end
	p5
end

# ╔═╡ Cell order:
# ╠═32e2bf36-e749-11f1-bcb7-a3a853c4441d
# ╠═eefe1cb6-90e3-4253-8640-017fe282c27c
# ╠═ea971183-d0c9-4148-9b31-2a09dd50e954
# ╠═539b8fd3-f82d-442f-8fb9-757757e31b3a
# ╠═b5eb1354-44c9-45a3-b9ad-3625a3498b6e
# ╠═7192fe57-1a4f-44cf-ae26-ecd123c20719
# ╠═15dfa605-2ae7-4c68-861d-c27810806929
# ╠═5e6157fa-c8f3-41e2-b760-cf1e93aae98a
# ╠═fdac42f3-0c2a-445b-b03c-7b07e3fb3730
# ╠═a95e67ab-a811-4c56-a2c8-2d72edf736e5
# ╠═f2d6c622-5e08-49f2-b5af-b3224ac57b21
# ╠═b6f47eef-d354-4cfb-811b-69ab91682894
# ╠═1a707304-9f9d-41ed-b6e0-a0e48636d67e
