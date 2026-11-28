### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ 35bb7e80-eb0c-11f1-9b55-b5b1cf728080
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
	using JpegTurbo
end

# ╔═╡ 7a1aaf6b-3d2a-49da-a106-44ec83ec1b2d
# Note: This code was adapted from https://discourse.julialang.org/t/get-dct-coefficients-of-jpeg-image/64640/5
begin
	# quantized DCT coefficients via libjpeg's jpeg_read_coefficients
	function read_dct(path)
		L = JpegTurbo.LibJpeg
		data = read(path)
		cinfo = Ref(L.jpeg_decompress_struct()); jerr = Ref{L.jpeg_error_mgr}()
		cinfo[].err = L.jpeg_std_error(jerr)
		L.jpeg_create_decompress(cinfo); L.jpeg_mem_src(cinfo, data, length(data))
		L.jpeg_read_header(cinfo, true)
		arrays = L.jpeg_read_coefficients(cinfo)
		access = unsafe_load(cinfo[].mem).access_virt_barray
		comps = map(1:cinfo[].num_components) do k
			ci = unsafe_load(cinfo[].comp_info, k)
			q = Int.(reshape(collect(unsafe_load(ci.quant_table).quantval), 8, 8)')
			coefs = zeros(Int16, 8, 8, ci.height_in_blocks, ci.width_in_blocks)
			for i in 1:ci.height_in_blocks
				row = unsafe_load(ccall(access, Ptr{Ptr{NTuple{64, Int16}}},
					(Ptr{Cvoid}, Ptr{Cvoid}, Cuint, Cuint, Cint),
					cinfo, unsafe_load(arrays, k), i - 1, 1, 0))
				for j in 1:ci.width_in_blocks
					coefs[:, :, i, j] = reshape(collect(unsafe_load(row, j)), 8, 8)'
				end
			end
			(q = q, coefs = coefs)
		end
		L.jpeg_destroy_decompress(cinfo)
		comps
	end

end

# ╔═╡ 977ea658-b0c4-4467-9585-8b2070b97710
begin
	valid_images = [1, 3, 4, 8, 9, 10, 12, 13, 15]
	coeffs = [(u, v) for u in 0:7 for v in 0:7]

	luma = Dict(n => read_dct(joinpath(@__DIR__, "rjpeg_plots",
		"grey$(lpad(n, 2, '0')).ppm.result.jpg"))[1] for n in valid_images)

	# Y[n, u, v] = DCT coefficient (u, v) of image n
	Y = Dict((n, u, v) => vec(luma[n].coefs[u+1, v+1, :, :])
		for n in valid_images for (u, v) in coeffs)

	# quantizations[n] = quality factor Q found by rjpeg
	quantizations = Dict(1 => 84, 3 => 58, 4 => 65, 8 => 80, 9 => 84,
		10 => 79, 12 => 76, 13 => 54, 15 => 85)
end

# ╔═╡ 7a1aaf6b-3d2a-49da-a106-44ec83ec1b2e
# histogram of the quantized (0,1) AC coefficient of image 1 across all blocks
histogram(Y[1, 0, 1], bins=-50.5:1:50.5, label="Y (0,1) coefficient")

# ╔═╡ ee9b26e6-383c-4660-846e-0aeaf4c7a3e3

function pf_spectrum(image, u, v; K = 30)
    k = Y[image, u, v]
    h = [count(==(n), k) for n in -K:K] ./ length(k)
    abs.(rfft(h))
end


# ╔═╡ dbd7d598-0e72-4cdd-900e-e3ce80f6670f
plot(pf_spectrum(1,0,1))

# ╔═╡ e0b78254-afe1-4ff9-961a-52a49a2cb6bc
let 
    ps = map(enumerate(valid_images)) do (i, img)
        p = plot(legend=false, title="Image $img")
        for (u,v) in coeffs
            (u, v) == (0, 0) && continue
            max(u, v) ≤ 2 || continue
            plot!(p, pf_spectrum(img,u,v))
        end
        p
    end
    plot(ps...; layout = (3, 3),
         plot_title = "PF spectrum")
end

# ╔═╡ Cell order:
# ╠═35bb7e80-eb0c-11f1-9b55-b5b1cf728080
# ╠═7a1aaf6b-3d2a-49da-a106-44ec83ec1b2d
# ╠═977ea658-b0c4-4467-9585-8b2070b97710
# ╠═7a1aaf6b-3d2a-49da-a106-44ec83ec1b2e
# ╠═ee9b26e6-383c-4660-846e-0aeaf4c7a3e3
# ╠═dbd7d598-0e72-4cdd-900e-e3ce80f6670f
# ╠═e0b78254-afe1-4ff9-961a-52a49a2cb6bc
