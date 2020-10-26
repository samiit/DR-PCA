### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ da5ad68e-1757-11eb-37ee-5ff74db1f6ca
using LinearAlgebra, Images, Plots

# ╔═╡ 6437c7fa-178e-11eb-3966-31c2f7ea3408
md"# Deconstructing PCA using a DR perspective

Results from the paper by [Narasimhan and Bhatt (2015)](https://www.sciencedirect.com/science/article/abs/pii/S0098135415000873)
"

# ╔═╡ 3a7ed398-1790-11eb-0662-07879f0beaf1
begin
	img = load("Flow_process.png")
	plot(img, axis=false)
end

# ╔═╡ 5661af4e-1791-11eb-0ded-3541c01717b3
md"### Constraint matrix"

# ╔═╡ 73ed6bc6-1759-11eb-024c-83cfae5b1b81
A = [1 1 -1 0 0 0;
	0 0 1 -1 0 0;
	0 0 0 1 -1 -1;
	0 -1 0 0 0 1]

# ╔═╡ 2fdcdcc0-175c-11eb-2a4e-5d913ceb7f52
md"## Measurements and Fluctuations
- Both, $F_1$ and $F_2$ are measured with fluctuations
- Base value for both is $10$
- Fluctuation are $\mathcal{N}(0,1)$ and $\mathcal{N}(0,2)$, respectively
- A total of $1000$ data points are available for the steady state flow balance equations
"

# ╔═╡ 1fb5fc10-175b-11eb-139c-8126aaca3a00
begin
	N = 1000
	F1_base, F1_fluct = 10, 1
	F2_base, F2_fluct = 10, 2
	F1ₜ = F1_base .+ (F1_fluct) * randn(N)
	F2ₜ = F2_base .+ (F2_fluct) * randn(N)
end;	

# ╔═╡ 7225218e-1791-11eb-2c02-617ad622b901
md"Breaking the constraint matrix into the measured $A_x$ and the unmeasured $A_u$ matrices."

# ╔═╡ d478db56-175e-11eb-31b6-d32eadc31c46
begin
	Aₓ = A[:,1:2]
	Aᵤ = A[:,3:end]
end;

# ╔═╡ bca0ca10-1791-11eb-3de6-b3fc37014809
md"For deriving the true values of the other variables, we need to substitute the true values of $F_1$ and $F_2$ (all $1000$ values) and obtain the corresponding values of the derived/dependent variables, $F_3$, $F_4$, $F_5$ and $F_6$.

$A_x [F_1, F_2]^T + A_u [F_3, F_4, F_5, F_6]^T = 0$
or

$A_u [F_3, F_4, F_5, F_6]^T = - A_x [F_1, F_2]^T$

"

# ╔═╡ ac79bc0a-175f-11eb-1205-df0aab8c88d6
b = -Aₓ*[F1ₜ F2ₜ]'

# ╔═╡ ea3993f8-175f-11eb-308a-11c6e6d08303
dependent_vars = Aᵤ\b

# ╔═╡ 7b461812-1760-11eb-1ca4-b9d591c18345
begin
	F3ₜ = dependent_vars[1,:]
	F4ₜ = dependent_vars[2, :]
	F5ₜ = dependent_vars[3, :]
	F6ₜ = dependent_vars[4, :]
end;

# ╔═╡ c06ad3ec-1760-11eb-178e-397e40a42dc4
size(F3ₜ), size(F4ₜ), size(F5ₜ), size(F6ₜ)

# ╔═╡ 78a1293e-1761-11eb-1d7c-67574000a6a4
md"#### Verify the true values"

# ╔═╡ 635e0f74-1761-11eb-0fc9-7b1c6bfc24df
X = [F1ₜ, F2ₜ, F3ₜ, F4ₜ, F5ₜ, F6ₜ]

# ╔═╡ 0590fbea-1761-11eb-0b0a-ad1cfe527582
A * X

# ╔═╡ 1f53b0c2-175c-11eb-0ad4-3b99740221e2
md"### SDE"

# ╔═╡ 089c4dde-175d-11eb-3702-3f3ad6fbc824
begin
	sde = [0.1, 0.08, 0.15, 0.2, 0.18, 0.1]
	Σ = diagm(sde.^2)
end

# ╔═╡ 9e3ed1e0-1762-11eb-3b40-0de3fe9e6c1f
md"#### Simulating measurements around the true values"

# ╔═╡ ac38bda6-1762-11eb-3c95-51a83eaa5e3d
Y = X + [sd * randn(N) for sd in sde]

# ╔═╡ 288f1b4e-175d-11eb-3c8b-1554fbd84e14
md"## Applying DR"

# ╔═╡ 45f72cbc-175d-11eb-0055-01c9815689c0
md"### LA route"

# ╔═╡ 5002007e-175d-11eb-2391-8fd69f92794a
W = I -Σ * A'* inv(A * Σ * A') * A

# ╔═╡ b44187ea-175e-11eb-2921-25e3c2345ca1
size(W)

# ╔═╡ c1e879ce-175e-11eb-0944-61396ca26a54
md"#### Estimates"

# ╔═╡ be181a9a-1761-11eb-1989-c1f15900ffdf
Xₑ = W * Y

# ╔═╡ 6ef7277e-176d-11eb-00ca-c95bc64fd8ed
md"##### Error covariance matrix"

# ╔═╡ 6a1c0092-1765-11eb-0ca7-7d6d5ee76b66
Σₑ = W * Σ * W'

# ╔═╡ 7d0befd4-176d-11eb-2881-69eea63cb4e9
diag(Σₑ)

# ╔═╡ 918e9662-1765-11eb-30b1-030cd8a03991
diag(Σₑ) .< diag(Σ)

# ╔═╡ 256ad438-1764-11eb-0783-874a2b0a8486
md"#### Verifying whether the estimates satisfy the constraints"

# ╔═╡ 19fb961c-1764-11eb-1eee-b9300189ab65
A * Xₑ

# ╔═╡ f09e3672-1764-11eb-3b4e-bff973f423b5
begin
	residual = reduce(hcat, A*Xₑ);
	size(residual)
end

# ╔═╡ 7ceeb166-1764-11eb-11ca-851e9d7fbdf6
all(x -> isapprox(x,0,atol=1e-12), residual)

# ╔═╡ 626ffc8a-1766-11eb-0830-4f2533596b99
md"#### Effect of DR"

# ╔═╡ 72f4527c-1766-11eb-3f89-03fd973eeaf3
begin
	rmse_before = [a.^2 for a in Y - X]
	rmse_before = [√sum(a) for a in rmse_before]
	
	rmse_after = [a.^2 for a in Xₑ - X]
	rmse_after = [√sum(a) for a in rmse_after]
end;

# ╔═╡ ba9ccd86-1767-11eb-0924-6154f763fa2e
rmse_before, rmse_after

# ╔═╡ a504e824-1768-11eb-0479-8b0fd218b00f
rmse_before .> rmse_after

# ╔═╡ 395f44de-1781-11eb-2f57-9326eb43c0c4
md"### PCA route"

# ╔═╡ 42e79e2a-1781-11eb-07c3-5595346dc7da
begin
	Y_mat = reduce(hcat, Y)'./ √N;
	size(Y_mat)
end

# ╔═╡ b7db36f8-1784-11eb-18e7-9df8eff09198
U, S, Vt = svd(Y_mat)

# ╔═╡ 35266c1c-1786-11eb-1dfc-89f1fdfc9f60
plot(S, marker="0", label="Singular values")

# ╔═╡ 946dd5a8-1785-11eb-3661-9734f2aaa832
threshold_var = 0.01;

# ╔═╡ 13c16910-1785-11eb-07b3-014ff0eb15bc
begin
	rel_var = 100*cumsum(reverse(S.^2))/sum(S.^2)
	plot(rel_var)
end

# ╔═╡ cb1d5724-1792-11eb-05ba-9bf4030ef2c0
rel_var

# ╔═╡ 66d7c00a-1787-11eb-2657-8f9e938a2346
m = findlast( rel_var .< threshold_var )

# ╔═╡ 2bc8d784-1788-11eb-3e11-99ee9641592b
begin
	n = size(A,2)
	p = n - m
end

# ╔═╡ c8986b5c-1787-11eb-2840-1d97394cdd34
md"### Denoised estimates
$X_e = √N (U_1  S_1  V_1 ^T)$
"

# ╔═╡ 96736630-1788-11eb-0528-c165152d4b6a
X_denoised = √N * U[:,1:p] * diagm(S[1:p]) * Vt[:,1:p]'

# ╔═╡ 15dec576-178a-11eb-1063-0fd1e64d75db
md"### Recovered Model

$A_m = {U_2}^T$
"

# ╔═╡ b68f2b3a-1789-11eb-18ed-c7a7668db14c
Aₘ = U[:,p+1:end]'

# ╔═╡ 88691d56-178a-11eb-178d-4b85c730b107
md"*Ofcourse, this estimated model is rotated w.r.t. to the true constraint matrix*

$A_m = QA$
But it is orthogonal to the measured estimates
"

# ╔═╡ 64ec59c2-1789-11eb-1b29-f7fe37cc910f
out_denoised = Aₘ*X_denoised

# ╔═╡ e1dc9cd4-178a-11eb-00c5-314ca4d343f7
all(x -> isapprox(x, 0, atol=1e-12), out_denoised)

# ╔═╡ d4299626-178c-11eb-2a26-d35ce8d672c9
md"**But does this satisfy the estimates obtained from LA based DR?**"

# ╔═╡ e1607cf6-178c-11eb-30c3-c50e55e2f919
out_der_model = Aₘ * Xₑ

# ╔═╡ ea6b8804-178c-11eb-310a-27b5f7374079
out_model = reduce(hcat, out_der_model)'

# ╔═╡ 0b6db504-178d-11eb-33cd-e13fdb6aa7dd
maximum(abs.(out_model)), minimum(abs.(out_model))

# ╔═╡ 90db70fa-178d-11eb-06b2-474370c4a531
100*maximum(abs.(out_model))/minimum(reduce(hcat, X)), 100*minimum(abs.(out_model))/minimum(reduce(hcat, X))

# ╔═╡ 5cd79d60-178d-11eb-0793-3329c24fcecb
md"*As is observed, the derived model from data alone, does not satisfy the estimates obtained from the traditional DR approach, but the residuals are still almost zero (upto 0.4%)*"

# ╔═╡ Cell order:
# ╟─6437c7fa-178e-11eb-3966-31c2f7ea3408
# ╠═da5ad68e-1757-11eb-37ee-5ff74db1f6ca
# ╟─3a7ed398-1790-11eb-0662-07879f0beaf1
# ╟─5661af4e-1791-11eb-0ded-3541c01717b3
# ╠═73ed6bc6-1759-11eb-024c-83cfae5b1b81
# ╟─2fdcdcc0-175c-11eb-2a4e-5d913ceb7f52
# ╠═1fb5fc10-175b-11eb-139c-8126aaca3a00
# ╟─7225218e-1791-11eb-2c02-617ad622b901
# ╠═d478db56-175e-11eb-31b6-d32eadc31c46
# ╟─bca0ca10-1791-11eb-3de6-b3fc37014809
# ╠═ac79bc0a-175f-11eb-1205-df0aab8c88d6
# ╠═ea3993f8-175f-11eb-308a-11c6e6d08303
# ╠═7b461812-1760-11eb-1ca4-b9d591c18345
# ╠═c06ad3ec-1760-11eb-178e-397e40a42dc4
# ╟─78a1293e-1761-11eb-1d7c-67574000a6a4
# ╠═635e0f74-1761-11eb-0fc9-7b1c6bfc24df
# ╠═0590fbea-1761-11eb-0b0a-ad1cfe527582
# ╟─1f53b0c2-175c-11eb-0ad4-3b99740221e2
# ╠═089c4dde-175d-11eb-3702-3f3ad6fbc824
# ╟─9e3ed1e0-1762-11eb-3b40-0de3fe9e6c1f
# ╠═ac38bda6-1762-11eb-3c95-51a83eaa5e3d
# ╟─288f1b4e-175d-11eb-3c8b-1554fbd84e14
# ╟─45f72cbc-175d-11eb-0055-01c9815689c0
# ╠═5002007e-175d-11eb-2391-8fd69f92794a
# ╠═b44187ea-175e-11eb-2921-25e3c2345ca1
# ╟─c1e879ce-175e-11eb-0944-61396ca26a54
# ╠═be181a9a-1761-11eb-1989-c1f15900ffdf
# ╟─6ef7277e-176d-11eb-00ca-c95bc64fd8ed
# ╠═6a1c0092-1765-11eb-0ca7-7d6d5ee76b66
# ╠═7d0befd4-176d-11eb-2881-69eea63cb4e9
# ╠═918e9662-1765-11eb-30b1-030cd8a03991
# ╟─256ad438-1764-11eb-0783-874a2b0a8486
# ╠═19fb961c-1764-11eb-1eee-b9300189ab65
# ╠═f09e3672-1764-11eb-3b4e-bff973f423b5
# ╠═7ceeb166-1764-11eb-11ca-851e9d7fbdf6
# ╟─626ffc8a-1766-11eb-0830-4f2533596b99
# ╠═72f4527c-1766-11eb-3f89-03fd973eeaf3
# ╠═ba9ccd86-1767-11eb-0924-6154f763fa2e
# ╠═a504e824-1768-11eb-0479-8b0fd218b00f
# ╟─395f44de-1781-11eb-2f57-9326eb43c0c4
# ╠═42e79e2a-1781-11eb-07c3-5595346dc7da
# ╠═b7db36f8-1784-11eb-18e7-9df8eff09198
# ╠═35266c1c-1786-11eb-1dfc-89f1fdfc9f60
# ╠═946dd5a8-1785-11eb-3661-9734f2aaa832
# ╠═13c16910-1785-11eb-07b3-014ff0eb15bc
# ╠═cb1d5724-1792-11eb-05ba-9bf4030ef2c0
# ╠═66d7c00a-1787-11eb-2657-8f9e938a2346
# ╠═2bc8d784-1788-11eb-3e11-99ee9641592b
# ╟─c8986b5c-1787-11eb-2840-1d97394cdd34
# ╠═96736630-1788-11eb-0528-c165152d4b6a
# ╟─15dec576-178a-11eb-1063-0fd1e64d75db
# ╠═b68f2b3a-1789-11eb-18ed-c7a7668db14c
# ╟─88691d56-178a-11eb-178d-4b85c730b107
# ╠═64ec59c2-1789-11eb-1b29-f7fe37cc910f
# ╠═e1dc9cd4-178a-11eb-00c5-314ca4d343f7
# ╟─d4299626-178c-11eb-2a26-d35ce8d672c9
# ╠═e1607cf6-178c-11eb-30c3-c50e55e2f919
# ╠═ea6b8804-178c-11eb-310a-27b5f7374079
# ╠═0b6db504-178d-11eb-33cd-e13fdb6aa7dd
# ╠═90db70fa-178d-11eb-06b2-474370c4a531
# ╟─5cd79d60-178d-11eb-0793-3329c24fcecb
