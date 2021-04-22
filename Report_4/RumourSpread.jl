### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ c26df480-a1a3-11eb-382c-95a5a19e6ff7
using GraphPlot, Plots, StatsBase, Random

# ╔═╡ dbe6e9f8-c977-47ef-9f46-d537c9de8e5b
html"<style>
main {
    max-width:75%;
    padding: unset;
    margin-right: unset !important;
    margin-bottom: 20px;
    align-self: center !important;
}
pluto-helpbox {
    visibility: hidden;
}
</style>"



# ╔═╡ 270fe308-9408-4ee3-82e8-8b1d871f76c6
import LightGraphs as lg

# ╔═╡ e9b9fdae-3ebc-454f-98e3-e6dee6efd27c
nC2(n) = n*(n-1)/2

# ╔═╡ 13b88f9c-c915-40c6-ab61-4da2b22d9b95
function mean_path_length(g)
	return sum(lg.floyd_warshall_shortest_paths(g).dists ) / (2 * nC2( length( lg.vertices(g) )))
end

# ╔═╡ 9110f100-2a1d-4659-a365-e664cc9139fb
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

# ╔═╡ 266e82bd-2c57-49d4-9161-1a45122a3c8b
x, y, z = plotCluster(1000, 10, vcat(0.001:0.001:0.01, 0.015:0.015:0.1, 0.115:0.1:1.0))

# ╔═╡ d751d227-9e09-42b4-b798-bfd54cefbc54
begin
	plot(x, y ./ y[1], color = :red, label = ["Clustering co-efficient (C/C0)" ""], xaxis = :log10, lw = 2,  st = [:line, :scatter])

	plot!(x, z ./ z[1], color = :blue, label = ["Mean path length (L/L0)" ""], xaxis = :log10, lw = 2, st = [:line, :scatter], legend = :best)	
end

# ╔═╡ f671b418-c77f-411e-adec-db13a1fe1c1d
g = lg.SimpleGraphs.watts_strogatz(10000, 10, 0.5)

# ╔═╡ 544187d9-8e8f-4380-9c19-9505edd22214
lg.degree_histogram(g)

# ╔═╡ aa89b9ef-cbe6-4115-8e43-db5c7654f762
plot(lg.degree_histogram(g), st = [:scatter, :line], lw = 2, color = :blue, label = ["degree distribution" ""])

# ╔═╡ 6caae0fc-a8c8-4781-acb1-b00286a6c031
function delayedStartModel(g, nNodes::Int64, nIter::Int64, spreadProb, delay::Int64, nAnti::Int64)
	
	infectedTime = zeros(Int64, nNodes) #number of steps spent in state 1
	rumourCount = zeros(Int64, nIter) #number of agents in state 1
	antiCount = zeros(Int64, nIter) #number of agents in state 2

	rumourState = zeros(Int64, nNodes)
	rumourState[rand(1:nNodes, 10)] .= 1 #initial rumour 
	oldState = copy(rumourState) #0: uninformed, 1: rumour, 2: anti-rumour

	for k in 1:nIter
		for vert in lg.vertices(g)
			if(oldState[vert] != 2)
				nbs = shuffle(lg.neighbors(g, vert))
				for nb in nbs
					rumourState[vert] = (rand() < spreadProb[oldState[nb] + 1]) ? oldState[nb] : oldState[vert]	
				end
				
				if(oldState[vert] == 1)
					infectedTime[vert] += 1
				end
				
			end
		end
		
		oldState = copy(rumourState)
		rumourCount[k] = count(x -> x==1, rumourState)
		antiCount[k] = count(x -> x==2, rumourState)

		if(k == delay)
			oldState[rand(1:nNodes, nAnti)] .= 2
		end
	end	
	
	return (rumourCount, antiCount, infectedTime)
end

# ╔═╡ 85b9b30b-de1a-4c25-823c-ef5ddea26cb1
begin
	nIter = 2000
	nNodes = 10000
	nEdges = 10
	p = 0.5
	gIn = lg.SimpleGraphs.watts_strogatz(nNodes, nEdges, p)
end

# ╔═╡ c507bc34-33c0-4dc4-a5e5-809856c445cc
begin
	probs = [0, 0.01, 0.03]
	delay = 500
	nAnti = 50
end

# ╔═╡ a8ac581c-8e42-4818-9226-83af84192d34
rC, aC, t = delayedStartModel(gIn, nNodes, nIter, probs, delay, nAnti)

# ╔═╡ 99712586-b968-4337-94f2-042a06bee02b
plot(1:nIter, [rC, aC, 10000 .- rC .- aC] , label = ["Rumour" "Anti-rumour" "Uninformed"], colour = [:red :blue :green], legend = :right, lw = 2)

# ╔═╡ 9a0ef478-b238-45b3-846f-da3c742efd0d
begin
	data1 = []
	data2 = []
	data3 = []
	for p in vcat(collect(0.01:0.02:0.1), collect(0.1:0.05:0.5))
		rC, aC, t = delayedStartModel(lg.SimpleGraphs.watts_strogatz(nNodes, nEdges, p), nNodes, nIter, probs, delay, nAnti)
		push!(data1, rC)
		push!(data2, aC)
		push!(data3, t)
	end
end

# ╔═╡ 350575cf-cda7-4e14-a0bc-64451e9f9525
begin
	plot(1:nIter, data1, lw=1, label = ["delay = 250" "delay = 500" "delay = 750" "delay = 1000" "delay = 1250" "delay = 1500"] ,legend = :topright)
	
	# xAx = collect(1:nIter)
	# maxData1 = findmax.(data1)
	# yMax = map(maxData1) do i; i[1]; end
	# xMax = map(maxData1) do i; xAx[i[2]]; end
	# plot!(xMax, yMax, lw = 2, color = :black, label = ["Maximum"])
end

# ╔═╡ c68c7200-ceb8-4acf-baab-37ddc83e5582
begin
	plot(1:nIter, data2, lw=1, label = ["delay = 250" "delay = 500" "delay = 750" "delay = 1000" "delay = 1250" "delay = 1500"] ,legend = :topleft)

end

# ╔═╡ 0d82b5d5-5636-47b4-9a27-2e6bd1db258a
begin
	tMax = [max(data3[i]...) for i in 1:length(data3)]
	plot(vcat(collect(0.01:0.02:0.1), collect(0.1:0.05:0.5)), tMax, lw=2, label = ["Max. time of rumour-belief vs. re-wiring probability" ""], legend = :topleft, st = [:line, :scatter], color = :blue)
end

# ╔═╡ c6cfbd45-de47-445f-ad61-7812524b59d3
begin
	tAv = [mean(data3[i]) for i in 1:length(data3)]
	plot(vcat(collect(0.01:0.02:0.1), collect(0.1:0.05:0.5)), tAv, lw=2, label = ["Avg. time of rumour-belief vs. re-wiring probability" ""], legend = :topleft, st = [:line, :scatter], color = :green, ylim = [0, 200])
end

# ╔═╡ 6da05dca-e144-4ece-b5e7-0725d5af1c06
function beaconModel(g, nNodes::Int64, nIter::Int64, spreadProb, nAnti::Int64)
	
	infectedTime = zeros(Int64, nNodes) #number of steps spent in state 1
	rumourCount = zeros(Int64, nIter) #number of agents in state 1
	antiCount = zeros(Int64, nIter) #number of agents in state 2

	rumourState = zeros(Int64, nNodes)
	rumourState[rand(1:nNodes, 10)] .= 1 #initial rumour 
	oldState = copy(rumourState) #0: uninformed, 1: rumour, 2: anti-rumour

	beacons = rand(1:nNodes, nAnti) #beacon indices
	
	for k in 1:nIter
		for vert in lg.vertices(g)
			if(oldState[vert] != 2)
				nbs = shuffle(lg.neighbors(g, vert))
				for nb in nbs
					rumourState[vert] = (rand() < spreadProb[oldState[nb] + 1]) ? oldState[nb] : oldState[vert]	
				end
				
				if(oldState[vert] == 1)
					infectedTime[vert] += 1
				end
				
				if((vert in beacons) && (1 in rumourState[lg.neighbors(g, vert)]))
					rumourState[vert] = 2
				end
			end
		end
		
		oldState = copy(rumourState)
		rumourCount[k] = count(x -> x==1, rumourState)
		antiCount[k] = count(x -> x==2, rumourState)
	end	
	
	return (rumourCount, antiCount, infectedTime)
end

# ╔═╡ a76b1fe9-7e41-4abd-ae74-ec3e2368a6ba
begin
	nIter1 = 2000
	nNodes1 = 10000
	nEdges1 = 10
	p1 = 0.5
	gIn1 = lg.SimpleGraphs.watts_strogatz(nNodes, nEdges, p)
end

# ╔═╡ 67da6b39-705b-42d0-b238-8f183c2c1a94
begin
	probs1 = [0, 0.01, 0.03]
	nAnti1 = 20
end

# ╔═╡ 1e521cfc-3964-4107-8340-a017fbc3fd81
a, b, c = beaconModel(gIn1, nNodes1, nIter1, probs1, nAnti1)

# ╔═╡ 1b6d1a56-bc4c-4d99-a5f2-da6f717bbcf1
plot(1:nIter1, [a, b, 10000 .- a .- b] , label = ["Rumour" "Anti-rumour" "Post" "Uninformed"], colour = [:red :blue :green], legend = :right)

# ╔═╡ 9e6adf8a-60bb-4801-8018-7d15721e4ff7
begin
	data11 = []
	data12 = []
	data13 = []
	for nAnti in 1:10
		rC, aC, t = beaconModel(gIn1, nNodes1, nIter1, probs1, nAnti)
		push!(data11, rC)
		push!(data12, aC)
		push!(data13, t)
	end
end

# ╔═╡ eae9e8b0-652d-4357-a22a-cfae382178c4
begin
	plot(1:nIter, data12, lw=1, label = ["#beacons = 1" "#beacons = 2" "#beacons = 3" "#beacons = 4" "#beacons = 5" "#beacons = 6" "#beacons = 7" "#beacons = 8" "#beacons = 9" "#beacons = 10"] ,legend = :topright)
	
	# xAx = collect(1:nIter)
	# maxData1 = findmax.(data1)
	# yMax = map(maxData1) do i; i[1]; end
	# xMax = map(maxData1) do i; xAx[i[2]]; end
	# plot!(xMax, yMax, lw = 2, color = :black, label = ["Maximum"])
end

# ╔═╡ 194a9c59-a031-4625-af7b-2b2f69fc9ad8
begin
	tMax1 = [max(data13[i]...) for i in 1:length(data13)]
	plot(1:10, tMax1, lw=2, label = ["Max. time of rumour-belief vs. #beacons" ""], legend = :topleft, st = [:line, :scatter], color = :blue)
end

# ╔═╡ 8e6dc070-b3ac-4d84-9db8-2b579ec579bb
begin
	tAv1 = [mean(data13[i]) for i in 1:length(data13)]
	plot(1:10, tAv1, lw=2, label = ["Avg. time of rumour-belief vs. #beacons" ""], legend = :topleft, st = [:line, :scatter], color = :green)
end

# ╔═╡ a84892e2-7e23-4509-a897-de148be7e968
function rumourModel(g, nNodes::Int64, nIter::Int64, spreadProb)
	
	rumourCount = zeros(Int64, nIter) #number of agents in state 1

	rumourState = zeros(Int64, nNodes)
	rumourState[rand(1:nNodes, 10)] .= 1 #initial rumour 
	oldState = copy(rumourState) #0: uninformed, 1: rumour, 2: anti-rumour

	for k in 1:nIter
		for vert in lg.vertices(g)
			if(oldState[vert] != 2)
				nbs = shuffle(lg.neighbors(g, vert))
				for nb in nbs
					rumourState[vert] = (rand() < spreadProb[oldState[nb] + 1]) ? oldState[nb] : oldState[vert]	
				end
			end
		end
		
		oldState = copy(rumourState)
		rumourCount[k] = count(x -> x==1, rumourState)
	end	
	
	return (rumourCount)
end

# ╔═╡ 7f32c561-109c-468d-8092-7f402e261e12
#gg = lg.SimpleGraphs.watts_strogatz(10000, 10, 0.5)
gg = lg.SimpleGraphs.watts_strogatz(10000, 10, 0.5)

# ╔═╡ 2ee734f4-24ae-48c0-b313-2051e69c2012
aa = rumourModel(gg, 10000, 2000, [0, 0.01, 0.03])

# ╔═╡ afc6fb22-4030-4756-9b75-4d5d375ae40a
plot(1:nIter1, [aa, 10000 .- aa] , label = ["Rumour" "Uninformed"], colour = [:red :blue :green], legend = :topright, lw=2)

# ╔═╡ 560bd9ae-b4bf-4aff-a0ae-e2fea736f984
begin
	plot(delete!(lg.degree_histogram(gg), 0), lw = 2, label = ["degree distribution" ""], st = [:line, :scatter])
end

# ╔═╡ Cell order:
# ╟─dbe6e9f8-c977-47ef-9f46-d537c9de8e5b
# ╠═c26df480-a1a3-11eb-382c-95a5a19e6ff7
# ╠═270fe308-9408-4ee3-82e8-8b1d871f76c6
# ╠═e9b9fdae-3ebc-454f-98e3-e6dee6efd27c
# ╠═13b88f9c-c915-40c6-ab61-4da2b22d9b95
# ╠═9110f100-2a1d-4659-a365-e664cc9139fb
# ╠═266e82bd-2c57-49d4-9161-1a45122a3c8b
# ╠═d751d227-9e09-42b4-b798-bfd54cefbc54
# ╠═f671b418-c77f-411e-adec-db13a1fe1c1d
# ╠═544187d9-8e8f-4380-9c19-9505edd22214
# ╠═aa89b9ef-cbe6-4115-8e43-db5c7654f762
# ╠═6caae0fc-a8c8-4781-acb1-b00286a6c031
# ╠═85b9b30b-de1a-4c25-823c-ef5ddea26cb1
# ╠═c507bc34-33c0-4dc4-a5e5-809856c445cc
# ╠═a8ac581c-8e42-4818-9226-83af84192d34
# ╠═99712586-b968-4337-94f2-042a06bee02b
# ╠═9a0ef478-b238-45b3-846f-da3c742efd0d
# ╠═350575cf-cda7-4e14-a0bc-64451e9f9525
# ╠═c68c7200-ceb8-4acf-baab-37ddc83e5582
# ╠═0d82b5d5-5636-47b4-9a27-2e6bd1db258a
# ╠═c6cfbd45-de47-445f-ad61-7812524b59d3
# ╠═6da05dca-e144-4ece-b5e7-0725d5af1c06
# ╠═a76b1fe9-7e41-4abd-ae74-ec3e2368a6ba
# ╠═67da6b39-705b-42d0-b238-8f183c2c1a94
# ╠═1e521cfc-3964-4107-8340-a017fbc3fd81
# ╠═1b6d1a56-bc4c-4d99-a5f2-da6f717bbcf1
# ╠═9e6adf8a-60bb-4801-8018-7d15721e4ff7
# ╠═eae9e8b0-652d-4357-a22a-cfae382178c4
# ╠═194a9c59-a031-4625-af7b-2b2f69fc9ad8
# ╠═8e6dc070-b3ac-4d84-9db8-2b579ec579bb
# ╠═a84892e2-7e23-4509-a897-de148be7e968
# ╠═7f32c561-109c-468d-8092-7f402e261e12
# ╠═2ee734f4-24ae-48c0-b313-2051e69c2012
# ╠═afc6fb22-4030-4756-9b75-4d5d375ae40a
# ╠═560bd9ae-b4bf-4aff-a0ae-e2fea736f984
