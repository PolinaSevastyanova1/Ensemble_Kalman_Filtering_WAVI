# using JLD2
# using EnsembleKalmanProcesses
# using NetCDF
# using Plots
# using Printf
# using Colors
# using Statistics

# @load "ensemble/output/eki.jld2" eki param_dict prior

# # initialising arrays
# iterations = length(eki.u)-1
# members = 20
# error_values = Array{Float64}(undef, iterations, members)

# #params stored in 6*10 array, of parameters and members, for each iteration
# weert, glen, bump, melt, pct, n = [Array{Float64}(undef, iterations, members) for _ in 1:6]

# # iterating over all iterations
# for i = 1:iterations
# 	u_iter = eki.u[i].stored_data
# 	weert[i, :], glen[i,:], bump[i,:], melt[i,:],pct[i,:], n[i,:] = eachrow(u_iter)
# end

# #param, plot title, ylabel, figure title, 
# params = [
# (weert, "Sliding Law Evolution", "weertman prefactor", "plots/weert_param_iters.png", (0.4,1.7)),
# (glen, "Glen's Flow Law Evolution", "glen prefactor", "plots/glen_param_iters.png", (0.4,1.7)),
# (bump, "Anthropogenic Bump Evolution", "bump prefactor", "plots/bump_param_iters.png", (-100,500)),
# (melt, "Melt Rate Prefactor Evolution", "melt rate prefactor", "plots/melt_param_iters.png", (0.5,15.5)),
# (pct, "Underlying Trend Evolution", "per century trend parameter", "plots/pct_param_iters.png", (-400,400)),
# (n, "n Exponent Evolution", "n exponent", "plots/n_exponent.png", (0.0,6.0))
# ]

# #and highlighting members 3,5,8,9 which led to the correct values. 
# #6 and 10 failed, so highlighting those as well, in dotted lines

# colours = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink,:red, :blue, :green, :orange, :purple, :cyan, :magenta, :yellow, :brown, :pink]
# for (param, title, ylabel, filename, ylims) in params

# 	plot(title=title, xlabel="Iteration no.", ylabel=ylabel, dpi=300, legend=:outerright, ylims=ylims)

# 	for member in 1:members
# 		plot!(1:iterations, param[:,member], label="Member $member", lw=1, color=colours[member])
# 	end

# 	savefig(filename)
# 	println("Plot saved in '$filename'")
# end




using JLD2
using EnsembleKalmanProcesses
using NetCDF
using Plots
using Statistics
using LaTeXStrings
using Printf
using Distributions
using StatsBase 
using Measures

# ──────────────────────────────────────────────────────────────────────────────
# 1) GLOBAL STYLING
# ──────────────────────────────────────────────────────────────────────────────
default(
  fontfamily   = "Computer Modern",
  titlefont    = font(14, "Computer Modern", :bold),
  guidefont    = font(14, "Computer Modern"),
  tickfont     = font(12, "Computer Modern"),
  legendfont   = font(12, "Computer Modern"),
  framestyle   = :box,
  grid         = false,
  dpi          = 1000,
)

# ──────────────────────────────────────────────────────────────────────────────
# 2) DATA FOLDERS & PARAMETERS
# ──────────────────────────────────────────────────────────────────────────────
const N_ITER = 20
const n_iter = 20
const N_MEMBERS = 20

# (symbol, title, ylabel, output filename, y‐limits)
params = [
  (:weert, "Sliding Law Evolution",        "Weertman prefactor",          "plots/n_weert.png", (0.4,1.7)),
  (:glen,  "Glen's Flow Law Evolution",    "Glen prefactor",              "plots/n_glen.png",  (0.4,1.7)),
  (:bump,  "Anthropogenic Bump Evolution", "Bump Amplitude (m)",              "plots/n_bump.png",  (-100,500)),
  (:melt,  "Melt Rate Prefactor Evolution","Melt rate prefactor",         "plots/n_melt.png",  (0.5,15.5)),
  (:pct,   "Underlying Trend Evolution",   "Anthropogenic trend (m/century)","plots/n_pct.png",(-300,300)),
  (:n_exponent,   "Glen Exponent Evolution",   "Glen 'n' exponent", "plots/n_exponent.png",(1.5,4.5))
]

# derive iteration count from first folder
@load "ensemble/output/eki.jld2" eki
obs = eki.observation_series.observations[1].samples


# ──────────────────────────────────────────────────────────────────────────────
# 3) PREPARE STORAGE FOR TRAJECTORIES & PRIORS
# ──────────────────────────────────────────────────────────────────────────────
all_failed  = Dict{Symbol,Array{Float64,2}}()
all_success = Dict{Symbol,Array{Float64,2}}()
all_prior   = Dict{Symbol,Vector{Float64}}()
for (sym,_,_,_,_) in params
  all_failed[sym]  = Array{Float64}(undef, N_ITER, 0)
  all_success[sym] = Array{Float64}(undef, N_ITER, 0)
  all_prior[sym]   = Float64[]
end

# ──────────────────────────────────────────────────────────────────────────────
# 4) ACCUMULATE ACROSS FOLDERS
# ──────────────────────────────────────────────────────────────────────────────


  @load "ensemble/output/eki.jld2" eki
  @load "ensemble/output/truth.jld2" y
  global obs = eki.observation_series.observations[1].samples
  μ_truth    = obs[1]
  tol        = 1.0
  obs_times = [285, 290, 295, 300]

  # ice‐volume to define success/fail
  ice_vals = Array{Float64}(undef, n_iter, N_MEMBERS)
  for i in 1:n_iter
    ice_vals[i,:] = vec(mean(eki.g[i].stored_data, dims=1))
  end
    # 1) grab the final‐iteration g matrix (obs × members)
    g_final = eki.g[end].stored_data  # size = (n_obs, N_MEMBERS)

    # 2) for each member, compare its entire column to y element‐wise
    #    and require all four to pass |pred – truth| ≤ 1
    mask_ok = [ all(abs.(g_final[:, m] .- y) .<= tol) for m in 1:N_MEMBERS ]

    # 3) indices
    ok_ix   = findall(mask_ok)
    fail_ix = findall(!, mask_ok)
  # unpack parameters
  weert = Array{Float64}(undef, n_iter, N_MEMBERS)
  glen  = Array{Float64}(undef, n_iter, N_MEMBERS)
  bump  = Array{Float64}(undef, n_iter, N_MEMBERS)
  melt  = Array{Float64}(undef, n_iter, N_MEMBERS)
  pct   = Array{Float64}(undef, n_iter, N_MEMBERS)
  n_exponent   = Array{Float64}(undef, n_iter, N_MEMBERS)
  for i in 1:n_iter
    row = eki.u[i].stored_data
    weert[i,:], glen[i,:], bump[i,:], melt[i,:], pct[i,:], n_exponent[i,:] = eachrow(row)
  end
  param_mats = Dict(:weert=>weert, :glen=>glen, :bump=>bump, :melt=>melt, :pct=>pct, :n_exponent=>n_exponent)

  #accumulate prior and final trajectories
  # for (sym, _, _, _, _) in params
  #   mat = param_mats[sym]
  #   append!(all_prior[sym], mat[1, :])                      # prior
  #   all_failed[sym]  = hcat(all_failed[sym],  mat[:, fail_ix])
  #   all_success[sym] = hcat(all_success[sym], mat[:, ok_ix])
  # end

    for (sym, _, _, _, _) in params
     mat = param_mats[sym]               # size = n_iter × N_MEMBERS

     # if fewer than N_ITER, pad by repeating the last row
     if n_iter < N_ITER
      lastrow = mat[end:end, :]         # 1×members
      extra   = repeat(lastrow, N_ITER - n_iter, 1)
      padded  = vcat(mat, extra)        # now N_ITER×members
     else
      padded  = mat
     end

     #prior always the iteration‐1 row
     append!(all_prior[sym], padded[1, :])

      # now safe to hcat the full N_ITER×? blocks
     all_failed[sym]  = hcat(all_failed[sym],  padded[:, fail_ix])
     all_success[sym] = hcat(all_success[sym], padded[:, ok_ix])
   end
# ──────────────────────────────────────────────────────────────────────────────
# 5) PLOT COMBINED: TRAJECTORIES + SHARED HISTOGRAM + GAUSSIANS SCALED TO COUNTS
# ──────────────────────────────────────────────────────────────────────────────
const FAIL_COLOR    = "#4A90E2"
const SUCCESS_COLOR = "#800020"
const NBINS = 20

for (sym, title, ylab, outfn, ylims) in params
  failed_mat  = all_failed[sym]
  success_mat = all_success[sym]
  prior_vals  = all_prior[sym]
  post_vals   = success_mat[end, :]

  # layout & share y-axis
  lyt = @layout([a{0.5w} b{0.5w}])
  p = plot(
  layout     = lyt,
  link       = :y,
  size       = (800,400),
  margin = 3mm  # add 5 mm all around
  )

  # (a) iteration trajectories
  plot!(
    p[1],
    1:N_ITER, failed_mat;
    color = FAIL_COLOR, alpha = 0.5, lw = 1, label = false
  )
  plot!(
    p[1],
    1:N_ITER, success_mat;
    color = SUCCESS_COLOR, alpha = 0.9, lw = 1, label = false
  )
  #title!(p[1], title)
  xlabel!(p[1], "Iteration")
  ylabel!(p[1], ylab)
  ylims!(p[1], ylims)

  clean_post = filter(!isnan, post_vals)
if length(clean_post)>1 && std(clean_post)>0
  h_post = fit(Normal, clean_post)
else
  @warn "No valid posterior samples for $sym – skipping gaussian overlay"
  h_post = nothing
end
  # (b) compute histogram counts
  edges = collect(range(ylims[1], ylims[2], length=NBINS+1))
  h_prior = fit(Histogram, prior_vals, edges; closed=:right)
  h_post  = fit(Histogram, post_vals,  edges; closed=:right)
  bw      = edges[2] - edges[1]

  # find the largest bin count
 max_count = maximum(vcat(h_prior.weights, h_post.weights))

 # pad it a bit (e.g. 5% headroom) and apply to the histogram panel
 xlims!(p[2], 0, max_count * 1.3)


  # prior histogram
  histogram!(
    p[2],
    prior_vals;
    orientation = :horizontal,
    bins        = edges,
    color       = FAIL_COLOR,
    alpha       = 0.15,
    label       = "Prior bins"
  )

#   # overlay fitted Gaussians scaled to counts
   xg = range(ylims[1], ylims[2], length=500)
   f_prior = fit(Normal, prior_vals)
   f_post  = fit(Normal, post_vals)
   scale_prior = length(prior_vals) * bw *1.8
   scale_post  = length(post_vals)  * bw *2.5#making the posterior bins twice as tall
#
   pdf_prior = pdf.(f_prior, xg) .* scale_prior
   pdf_post  = pdf.(f_post,  xg) .* scale_post

  #Fitted prior Guassians
  plot!(
    p[2],
    pdf_prior, xg;
    lw    = 2,
    color = FAIL_COLOR,
    label = false,
    alpha = 1
  )
#  scale_post = length(post_vals) * bw * 2.5       # bump Gaussian 20% larger
  pdf_post   = pdf.(f_post, xg) .* scale_post

#  #Fitted posterior gaussian
 plot!(
  p[2],
  pdf_post, xg;
  lw    = 2,
  color = SUCCESS_COLOR,
  label = false,
  alpha = 1
 )

#  # compute ±1σ limits
  μ_prior, σ_prior = mean(prior_vals), std(prior_vals)
  μ_post,  σ_post  = mean(post_vals),  std(post_vals)

#  # compute PDF heights at those y‐positions
  pdf_at(μ, f, scale) = pdf(f, μ) * scale

  h_prior = fit(Normal, prior_vals)
  h_post  = fit(Normal, post_vals)

 hline_vals = [
  (μ_prior,        FAIL_COLOR),
  (μ_prior+σ_prior,FAIL_COLOR),
  (μ_prior-σ_prior,FAIL_COLOR),
  (μ_post,         SUCCESS_COLOR),
  (μ_post+σ_post,  SUCCESS_COLOR),
  (μ_post-σ_post,  SUCCESS_COLOR)
 ]

# — Prior ±1σ lines
for yv in (μ_prior, μ_prior-σ_prior, μ_prior+σ_prior)
    sc = length(prior_vals) * bw
    xh = pdf(f_prior, yv) * sc
    plot!(
      p[2],
      [0, xh*1.8], [yv, yv];
      color = FAIL_COLOR,
      lw    = 2,
      label = false
    )
end

 # mask indices within ±σ
 mask_p = (xg .>= μ_prior - σ_prior) .& (xg .<= μ_prior + σ_prior)
 mask_q = (xg .>= μ_post  - σ_post ) .& (xg .<= μ_post  + σ_post)

 # build and plot a filled shape for the prior
 ys_p  = xg[mask_p]
 xs_p  = pdf_prior[mask_p]
 poly_x = [0.0; xs_p; 0.0]
 poly_y = [μ_prior-σ_prior; ys_p; μ_prior+σ_prior]
 plot!(
  p[2],
  poly_x, poly_y;
  seriestype = :shape,
  fillcolor  = FAIL_COLOR,
  fillalpha  = 0.15,
  linealpha  = 0,
  label      = false
 )


#build and plot a filled shape for the posterior
 mask_q  = (xg .>= μ_post-σ_post) .& (xg .<= μ_post+σ_post)
 ys_q   = xg[mask_q]
 xs_q   = pdf_post[mask_q]
 poly_x2 = [0.0; xs_q; 0.0]
 poly_y2 = [μ_post-σ_post; ys_q; μ_post+σ_post]
 plot!(
  p[2],
  poly_x2, poly_y2;
  seriestype = :shape,
  fillcolor  = "#d46c7b",
  fillalpha  = 0.8,
  linealpha  = 0,
  label      = false
 )


  # posterior histogram
 histogram!(
  p[2],
  post_vals;
  orientation = :horizontal,
  bins        = edges,
  weights     = fill(1.0, length(post_vals)),       # scale counts ×2
  color       = SUCCESS_COLOR,                          # lighter burgundy
  alpha       = 0.9,                                # opaque
  label       = "Posterior bins"
 )

 #— Posterior ±1σ lines
for yv in (μ_post, μ_post-σ_post, μ_post+σ_post)
    sc = length(post_vals) * bw
    xh = pdf(f_post, yv) * sc
    plot!(
      p[2],
      [0, xh*2.5], [yv, yv];
      color = SUCCESS_COLOR,
      lw    = 2,
      label = false
    )
end

 # 2) prepare info string (use LaTeXStrings if you like)
 info = @sprintf(
  "Prior: μ=%.2f, σ=%.2f\nPosterior: μ=%.2f, σ=%.2f",
  μ_prior, σ_prior, μ_post, σ_post
 )

 # 3) decide a data‐coordinate corner in panel p[2]:
 #    x at 90% of horizontal span, y at 10% up from y_min
 xpos = max_count * 0.9
 ypos = ylims[1] + 0.1*(ylims[2] - ylims[1])

 #dummy prior entry
 plot!(
  p[2],
  [NaN], [NaN];
  color = FAIL_COLOR,
  lw    = 2,
  label = @sprintf("Prior fit: μ=%.2f, σ=%.2f", μ_prior, σ_prior)
 )

#dummy posterior entry
 plot!(
  p[2],
  [NaN], [NaN];
  color = SUCCESS_COLOR,
  lw    = 2,
  label = @sprintf("Posterior fit: μ=%.2f, σ=%.2f", μ_post, σ_post)
 )

 plot!(
  p[2],
  legend            = :bottomright,
  legendfontsize    = 8,
  legendframe       = :box,
  legendbackgroundcolor = :white,
  legendbordercolor = :black,
  legendframealpha  = 0.9
 )
  xlabel!(p[2], "Count")
  ylabel!(p[2], ylab)
  ylims!(p[2], ylims)

  savefig(p, outfn)
  println("Saved combined figure for $sym → $outfn")
end
