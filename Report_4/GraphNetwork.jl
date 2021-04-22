### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 444a0796-a101-11eb-2f1c-77402a2f39f4
using GraphPlot, Plots, StatsBase

# ╔═╡ 2795feb6-59e8-4f4e-85a2-9ec45b9163b4
import LightGraphs as lg

# ╔═╡ bd5e27eb-16be-4d26-806d-faa028f00047
nC2(n) = n*(n-1)/2

# ╔═╡ 6ce026ec-7c00-4c20-90cf-cc6d6fbe3d7a
function mean_path_length(g)
	return sum(lg.floyd_warshall_shortest_paths(g).dists ) / (2 * nC2( length( lg.vertices(g) )))
end

# ╔═╡ cb7b1b1c-698c-4767-ae9b-1524de594d84
function plotCluster(nodes, Nedges, urange)
	nSteps = length(urange)
	clustering = zeros(Float64, nSteps)
	mpl = zeros(Float64, nSteps)
	
	for (i, p) in enumerate(urange)
		
		g = lg.SimpleGraphs.watts_strogatz(nodes, Nedges, p)
		
		clustering[i] = lg.global_clustering_coefficient(g)
		mpl[i] = mean_path_length(g)
	end
	
	return (urange, clustering, mpl)
end

# ╔═╡ bae8c338-70bc-4cb4-9e78-beec1533a3c9
x, y, z = plotCluster(1000, 10, vcat(0.001:0.001:0.01, 0.015:0.015:0.1, 0.115:0.1:1.0))

# ╔═╡ f0f99c3f-e76b-45bb-b24f-b7f7fa45c2bf
begin
	plot(x, y ./ y[1], color = :red, label = ["Clustering co-efficient" ""], xaxis = :log10, lw = 2,  st = [:line, :scatter])

	plot!(x, z ./ z[1], color = :blue, label = ["Mean path length" ""], xaxis = :log10, lw = 2, st = [:line, :scatter], legend = :bottomleft)	
end

# ╔═╡ 94e2e470-0dca-40c4-9fb8-05ae48cdedeb
g = lg.SimpleGraphs.barabasi_albert(5000, 100, 20)

# ╔═╡ a19637c4-822d-4fec-8f76-6bf895d0efca
begin
	n = lg.degree_histogram(g)
	delete!(n, 0)
end

# ╔═╡ 3b6eb5a1-c167-4a6c-b699-7e5f9a013b68
plot(n, xaxis=:log10, yaxis =:log10, st = [:scatter])

# ╔═╡ ecc4a252-955b-49d1-badc-cf6970f0d325
function competitionModel(g, nNodes, r, nIter)
	s = rand(0:1, nNodes) #strategies of agents
	sOld = copy(s)
	payoff = zeros(Float64, nNodes) #payoff for each agent after a set of games
	steadyState = nIter ÷ 100
	coop_density = zeros(Float64, steadyState)
	j = 1 #loop counter for density	
	
	for i ∈ 1:nIter
		payoff .= 0 #reset payoffs
		
		for vert in lg.vertices(g)
			nbs = vcat(vert, lg.neighbors(g, vert))
			payoff[nbs] .-= s[nbs]
			payoff[nbs] .+= (r/(length(nbs)) * sum(s[nbs]))
		end
		
		for vert in lg.vertices(g)
			random_nb = rand(lg.neighbors(g, vert))
			
			if(payoff[vert] < payoff[random_nb])
				s[vert] = sOld[random_nb]				
			end
		end
		
		if(i > nIter - steadyState)
			coop_density[j] = count(x -> x==1, s)
	 		j += 1
		end
		
		sOld = copy(s)
	end
	return mean(coop_density ./ nNodes)
end

# ╔═╡ b79fd39c-8811-4dd5-9f5e-bb64394c0a9f
begin
	nNodes = 1000
	nEdges = 6
	p = 0.05
	r = 5
	nIter = Int(1e5)
	
	gSnW = lg.SimpleGraphs.watts_strogatz(nNodes, nEdges, p) #network
	
	(lg.global_clustering_coefficient(gSnW), mean_path_length(gSnW))
end

# ╔═╡ 3beafb96-e26f-4e07-a8da-03d12f5b84b3
begin
	rVals = 5.5:0.25:7.5
	densities = zeros(length(rVals)) 
	
	for (i, r) in enumerate(rVals)
		densities[i] = competitionModel(gSnW, nNodes, r, nIter)
	end
end

# ╔═╡ 0a1561b4-cbe6-48d9-b781-f98d89b10f17
plot(rVals, densities, st = [:line, :scatter], color = :blue, ylims = [0, 1])

# ╔═╡ a406a296-0fa2-4786-a5b9-0f40b0a216c2
function competitionModelData(g, nNodes, r, nIter)
	s = rand(0:1, nNodes) #strategies of agents
	sOld = copy(s)
	payoff = zeros(Float64, nNodes) #payoff for each agent after a set of games
	steadyState = nIter ÷ 100
	coop_density = zeros(Float64, nIter)
	maxP = 0.0
	
	for i ∈ 1:nIter
		payoff .= 0 #reset payoffs
		coop_density[i] = count(x -> x==1, s)
		
		for vert in lg.vertices(g)
			nbs = vcat(vert, lg.neighbors(g, vert))
			payoff[nbs] .-= s[nbs]
			payoff[nbs] .+= (r/(length(nbs)) * sum(s[nbs]))
		end
		
		maxP = max(maxP, payoff...)
		
		for vert in lg.vertices(g)
			random_nb = rand(lg.neighbors(g, vert))
			
			if(rand() < (payoff[vert] - payoff[random_nb])/maxP)
				s[vert] = sOld[random_nb]				
			end
		end
		
		sOld = copy(s)
	end
	return coop_density ./ nNodes
end

# ╔═╡ 3e55a0d8-745c-484d-9966-8c5afeed6b64
begin
	nNodes1 = 1000
	nEdges1 = 6
	p1 = 0.1
	r1 = 5
	nIter1 = Int(1e5)
	
	gSnW1 = lg.SimpleGraphs.watts_strogatz(nNodes1, nEdges1, p1) #network
	
	(lg.global_clustering_coefficient(gSnW1), mean_path_length(gSnW1))
end

# ╔═╡ e7efbc4d-5e2b-423f-b96d-e3a042d7585e
dens = competitionModelData(gSnW1, nNodes1, r1, nIter1)

# ╔═╡ bc62662e-bad2-4361-ac91-a5c731b01135
plot(1:length(dens), dens, ylim = [0,1], xaxis = :log10)

# ╔═╡ 46cd0387-5d70-416a-9927-3c05440dd7fd
max.([1, 3, 5]...)

# ╔═╡ Cell order:
# ╠═444a0796-a101-11eb-2f1c-77402a2f39f4
# ╠═2795feb6-59e8-4f4e-85a2-9ec45b9163b4
# ╠═bd5e27eb-16be-4d26-806d-faa028f00047
# ╠═6ce026ec-7c00-4c20-90cf-cc6d6fbe3d7a
# ╠═cb7b1b1c-698c-4767-ae9b-1524de594d84
# ╠═bae8c338-70bc-4cb4-9e78-beec1533a3c9
# ╠═f0f99c3f-e76b-45bb-b24f-b7f7fa45c2bf
# ╠═94e2e470-0dca-40c4-9fb8-05ae48cdedeb
# ╠═a19637c4-822d-4fec-8f76-6bf895d0efca
# ╠═3b6eb5a1-c167-4a6c-b699-7e5f9a013b68
# ╠═ecc4a252-955b-49d1-badc-cf6970f0d325
# ╠═b79fd39c-8811-4dd5-9f5e-bb64394c0a9f
# ╠═3beafb96-e26f-4e07-a8da-03d12f5b84b3
# ╠═0a1561b4-cbe6-48d9-b781-f98d89b10f17
# ╠═a406a296-0fa2-4786-a5b9-0f40b0a216c2
# ╠═3e55a0d8-745c-484d-9966-8c5afeed6b64
# ╠═e7efbc4d-5e2b-423f-b96d-e3a042d7585e
# ╠═bc62662e-bad2-4361-ac91-a5c731b01135
# ╠═46cd0387-5d70-416a-9927-3c05440dd7fd
