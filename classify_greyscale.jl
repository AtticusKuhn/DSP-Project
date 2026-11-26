### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ cdeef0a5-bbbb-4146-9de9-66f09e6e8d1c
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


# ╔═╡ 2511b0aa-e7f7-11f1-8570-1da92e9b9b0a
images = [Float64.(rawview(channelview(load("./dataset/grey$(lpad(q, 2, '0')).png")))) for q in 1:15]

# ╔═╡ 5f4cb233-c477-4558-902f-a990f0137606
Z(x) = abs.(x[1:end-1, 1:end-1] .- x[2:end, 1:end-1] .- x[1:end-1, 2:end] .+ x[2:end, 2:end])

# ╔═╡ 41ee7ebb-ee81-4459-b032-7f4ea18d9a5c
Z1_Z2(image) = let 
	Z_png = Z(image)
	H, W = size(Z_png)
	p = mod1.(1:H, 8)            
	q = mod1.(1:W, 8)            
	
	Z1 = Z_png[p .== 4, q .== 4]      
	Z2 = Z_png[p .== 8, q .== 8]   
	(Z1, Z2)
end;

# ╔═╡ dc8b2699-e0ad-488a-bd77-6dd7ac32ea64
# size(Z(png))

# ╔═╡ 7b4c6a68-89e5-4f18-b158-f1e1e5cbab3d
    zhist(z) = [count(==(v), z) for v in 0:510] ./ length(z)

# ╔═╡ 5edb9511-5e47-43d9-afcf-fcb8ec2f066c
  d(img) = (z = Z1_Z2(img); cumsum(abs.(zhist(z[1]) .- zhist(z[2]))))

# ╔═╡ b71150ce-b273-4ede-809a-f109ddc2d215
let
        p = plot(xlabel="Z", ylabel="|H_I − H_II|")
      

        for (i, img) in enumerate(images)
                plot!(p, 0:510, d(img), label="image $i")
        end
        p
end


# ╔═╡ cc6ab6ed-6cb0-46e6-9299-36f2d5f3094c
Zscores = [(z = Z1_Z2(img); sum(abs.(zhist(z[1]) .- zhist(z[2])))) for img in images]

# ╔═╡ be3f18b9-9816-42ba-9412-4779baaabf6b
	  c_hist(c, vals) = [count(==(v), round.(Int, c)) for v in vals] ./ length(c)


# ╔═╡ a08ad2d8-4845-495f-92b5-b650486b84dc
begin
	blocks(x) = [x[i:i+7, j:j+7] for i in 1:8:8(size(x,1)÷8), j in 1:8:8(size(x,2)÷8)];
	Y(image) = [dct(B) for B in blocks(image)]
	cuv(image, u, v) = [y[u+1, v+1] for y in Y(image)]
	c01(image) = [y[1, 2] for y in Y(image)]
end

# ╔═╡ d96d1ae6-a6e6-43cb-92aa-41ad4f291069

function pf_spectrum(image, u, v, q2; K = 30)
    k = round.(Int, cuv(image, u, v) ./ q2)
    h = [count(==(n), k) for n in -K:K] ./ length(k)
    abs.(rfft(h))
end


# ╔═╡ 52bff5df-9a9d-4bce-85ba-11fd7e9b9332
begin
	vals = -60:60
	p3 = plot(label="uncompressed", title="Histogram of Y(0,1)")
	for (i, img) in enumerate(images) 
		plot!(p3, abs.(rfft(c_hist(c01(img), vals))), label="image $i")
	end
	p3
end

# ╔═╡ 3b96999b-a311-400e-be80-efaf095fa536
begin
	r(image,s ,m, n, q) = round(Y(image)[m,n][s] / q)
	# L = log likelhood
	D(m) = m == 0 || m == 4 ? 2 : m == 2 || m == 6 ? 2cos(pi/4) : 2cos(pi/8)cos(pi/8)
	B(m, n) = D(m)D(n)
	G(x) = exp( - 6 * x^2)
	function L(c, m, n, q)                     
        b = B(m, n)
        K = div(b, q)
        sum(log(sum(G(y - round(y/q)*q - k*q) for k in -K:K))  # Σ_s Σ_k
            for y in c) + length(c) * log(q)                                    # + N log q
	end
end


# ╔═╡ d4364b76-a452-4f31-a4aa-cb4f719849b0
begin
	block_is_good(b) = maximum(b) < 255 && minimum(b) > 0 && maximum(b) > minimum(b)
	  coeffs(image, m, n) = [dct(blk .- 128)[m+1, n+1] for blk in blocks(image) if block_is_good(blk)]
	qhat(image, m, n) = (c = coeffs(image, m, n); argmax(q -> L(c, m, n, q), 1:60))
end

# ╔═╡ 5b45d4ca-05ea-4adf-9e6b-58bc7f1c1048
begin
	p6 = plot(label="q=q", title="Histogram of Y(0,1)")
	for (i, img) in enumerate(images)
		q_esimate = qhat(img, 0, 1)
		plot!(p6, pf_spectrum(img, 0, 1, q_esimate), label="image $i")
	end
	p6
end

# ╔═╡ b525f8e1-055a-4db9-8474-12b2a3542609
q_estimates = qhat.(images, 0, 1)


# ╔═╡ e719cd99-09e5-489c-b18d-ba0272c566cb
begin
        # IJG standard luminance quantization table
        const T_LUM = [16  11  10  16  24  40  51  61
                       12  12  14  19  26  58  60  55
                       14  13  16  24  40  57  69  56
                       14  17  22  29  51  87  80  62
                       18  22  37  56  68 109 103  77
                       24  35  55  64  81 104 113  92
                       49  64  78  87 103 121 120 101
                       72  92  95  98 112 100 103  99]

        # The table cjpeg -quality Q uses
        function qtable(Q)
                S = Q < 50 ? 5000 ÷ Q : 200 - 2Q
                clamp.((T_LUM .* S .+ 50) .÷ 100, 1, 255)
        end

        Q_candidates(q, m, n) = [Q for Q in 1:100 if qtable(Q)[m+1, n+1] == q]

        # Several estimates at once
        Q_candidates(ests) =
                reduce(intersect, (Q_candidates(q, m, n) for ((m, n), q) in ests); init = 1:100)
end

# ╔═╡ 250e311a-4299-4e91-a73b-a8039dc8cfe3
begin
	
	coefs = [(0,0), (0,1), (1,0), (1,1), (0,2), (2,0)]
	Q_est = [Q_candidates(Dict(c => qhat(img, c...) for c in coefs)) for img in images]
	
end

# ╔═╡ 752f769a-301e-4777-b6c4-2d6e1aff3c8e
show(Q_est)

# ╔═╡ 4e2d8c6f-ecc0-46bf-b677-a455078683cd

begin
        alpha_to_Q(α) = α <= 1 ? 100 - 50α : 50 / α

        function Qstar(img; coefs = [(0,0), (0,1), (1,0), (1,1), (0,2), (2,0), (1,2), (2,1), (2,2)])
                qs = [qhat(img, m, n) for (m, n) in coefs]
                ts = [T_LUM[m+1, n+1] for (m, n) in coefs]
                # rel = qs .> 1
                # any(rel) || return (Q = missing, per_coef = Float64[], spread = 0.0)  # all q̂ = 1 → never compressed
                perQ = alpha_to_Q.(qs ./ ts)           
                α̂ = sum(qs) / sum(ts)                  
                (Q = alpha_to_Q(α̂), per_coef = perQ, spread = maximum(perQ) - minimum(perQ))
        end
end

# ╔═╡ 50bbc063-affe-4cbd-8cfb-d5c24e614d07
qstars = Qstar.(images)

# ╔═╡ cb2a7765-e1a2-45f7-87db-ea55bddf3df9
show(qstars)

# ╔═╡ 4fa87e4e-4f47-403a-a608-113a611b8b53
begin
	p7 = plot(label="q=q", title="Histogram of various coefs")
	for c in coefs
		for (i, img) in enumerate(images)
			q_esimate = qhat(img, c[1], c[2])
			plot!(p7, pf_spectrum(img, c[1], c[2], q_esimate), label="image $i, c=$c", color=i)
		end
	end
	p7
end

# ╔═╡ dcb2ec4d-b19d-44ec-9318-c9838dc17a7c
begin
    function peak_score(img, m, n; K = 30, w = 3, cutoff = 5)
        q = qhat(img, m, n)
        s = pf_spectrum(img, m, n, q; K)
        # prominence over a symmetric neighbourhood; cutoff only restricts which k are candidates
        prom = [s[k] - median(s[k-w:min(end, k+w)]) for k in cutoff:length(s)]
        k = argmax(prom)
        (q = q, idx = k + cutoff - 1, score = prom[k])
    end
    
    ac_coefs = [(0,1), (1,0), (1,1), (0,2), (2,0), (1,2), (2,1), (2,2)]
    peaks = Dict((i, c) => peak_score(img, c...) for (i, img) in enumerate(images) for c in ac_coefs)
    
end

# ╔═╡ ae47021d-4711-4a96-a4d2-068a277d97da
peaks

# ╔═╡ 78c8609e-41f8-4714-b7b6-ef5e576bb73d
show(peaks)

# ╔═╡ 7fa006c4-f10d-4d19-94d2-9d6ac0b6826c
let coef = (0, 1), K = 15
    m, n = coef
    ps = map(enumerate(images)) do (i, img)
        q = qhat(img, m, n)
        k = round.(Int, cuv(img, m, n) ./ q)
        h = [count(==(j), k) for j in -K:K] ./ length(k)
        bar(-K:K, h; legend = false, linewidth = 0,
            title = "grey$(lpad(i, 2, '0'))  q2=$q  score=$(round(peaks[(i, coef)].score, digits = 4))",
            titlefontsize = 8, tickfontsize = 6)
    end
    plot(ps...; layout = (5, 3), size = (1100, 1300),
         plot_title = "Histogram of round(c$(m)$(n) / q2)")
end

# ╔═╡ b8bc5f1c-e8b5-47cb-90de-356c63b577ef
let K = 15
    ps = map(enumerate(images)) do (i, img)
        p = plot(title = "grey$(lpad(i, 2, '0'))", yscale = :log10, titlefontsize = 9, tickfontsize = 6, legend = false)
        for c in ac_coefs
            q = peaks[(i, c)].q       
            k = round.(Int, cuv(img, c...) ./ q)
            h = [count(==(j), k) for j in -K:K] ./ length(k)
            plot!(p, -K:K, h; label = "$c q=$q")
        end
        p
    end
    plot(ps...; layout = (5, 3), size = (1100, 1300),
         plot_title = "Histograms of round(c_uv / q2), all AC coefficients")
end

# ╔═╡ Cell order:
# ╠═cdeef0a5-bbbb-4146-9de9-66f09e6e8d1c
# ╠═2511b0aa-e7f7-11f1-8570-1da92e9b9b0a
# ╠═5f4cb233-c477-4558-902f-a990f0137606
# ╠═41ee7ebb-ee81-4459-b032-7f4ea18d9a5c
# ╠═dc8b2699-e0ad-488a-bd77-6dd7ac32ea64
# ╠═7b4c6a68-89e5-4f18-b158-f1e1e5cbab3d
# ╠═5edb9511-5e47-43d9-afcf-fcb8ec2f066c
# ╠═b71150ce-b273-4ede-809a-f109ddc2d215
# ╠═cc6ab6ed-6cb0-46e6-9299-36f2d5f3094c
# ╠═be3f18b9-9816-42ba-9412-4779baaabf6b
# ╠═d96d1ae6-a6e6-43cb-92aa-41ad4f291069
# ╠═a08ad2d8-4845-495f-92b5-b650486b84dc
# ╠═52bff5df-9a9d-4bce-85ba-11fd7e9b9332
# ╠═3b96999b-a311-400e-be80-efaf095fa536
# ╠═d4364b76-a452-4f31-a4aa-cb4f719849b0
# ╠═5b45d4ca-05ea-4adf-9e6b-58bc7f1c1048
# ╠═b525f8e1-055a-4db9-8474-12b2a3542609
# ╠═e719cd99-09e5-489c-b18d-ba0272c566cb
# ╠═250e311a-4299-4e91-a73b-a8039dc8cfe3
# ╠═752f769a-301e-4777-b6c4-2d6e1aff3c8e
# ╠═4e2d8c6f-ecc0-46bf-b677-a455078683cd
# ╠═50bbc063-affe-4cbd-8cfb-d5c24e614d07
# ╠═cb2a7765-e1a2-45f7-87db-ea55bddf3df9
# ╠═4fa87e4e-4f47-403a-a608-113a611b8b53
# ╠═dcb2ec4d-b19d-44ec-9318-c9838dc17a7c
# ╠═ae47021d-4711-4a96-a4d2-068a277d97da
# ╠═78c8609e-41f8-4714-b7b6-ef5e576bb73d
# ╠═7fa006c4-f10d-4d19-94d2-9d6ac0b6826c
# ╠═b8bc5f1c-e8b5-47cb-90de-356c63b577ef
