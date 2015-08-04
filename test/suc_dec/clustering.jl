# Graph clustering

using PyPlot;
using Graphs;
using GraphViz;
using MetisMod;

COLORS = ["blue", "red", "green", "orange", "purple", "yellow", "indigo", "pink", "aquamarine", "brown", "coral", "cyan", "darkkhaki", "gray"];

# Construct a graph based on the network and partition it
function build_and_cluster_graph(nBuses, flowin_bus_id, flowout_bus_id, weights=Dict{(Int,Int), Int}(), nclusters=3)
	graph, edgelist = graph_from_data(nBuses, flowin_bus_id, flowout_bus_id);

	# Unweighted version if weight dict is empty
	if isempty(weights)
		weights = Dict{(Int,Int), Int}();
		for edge in edgelist
			weights[edge] = 1;
		end
	end

	partition = cluster(graph, weights, nclusters, :metis_kway);

	# TODO: Temporary, for testing
	# for i in 1:length(partition)
	# 	partition[i] = 1
	# end
	# partition[10] = 2

	plot_partitioned_graph(nBuses, edgelist, partition);

	# Print info
	counts = zeros(Int, nclusters);
	for p in partition
		counts[p] += 1;
	end
	print("Partition sizes: ");
	println(counts);
	print("Partition: ");
	println(partition);

	return graph, edgelist, partition
end

# Collect cut edges from a partition (i.e. edges with endpoints in different classes)
function collect_cut_edges(nBuses, edgelist, partition)
	cut_edges = Int[];
	for (i, edge) in enumerate(edgelist)
		assert(edge[1] != edge[2]);
		if partition[edge[1]] != partition[edge[2]]
			push!(cut_edges, i);
		end
	end
	return cut_edges;
end

# Plot graph with arguments to neato
function plot_args(g::AbstractGraph, opts::String="")
	if isempty(opts)
		stdin, proc = open(`neato -Tx11`, "w")
	else
		stdin, proc = open(`neato -Tx11 $opts`, "w")
	end
    d = AttributeDict()
	to_dot(g, stdin, d);
	close(stdin)
end

# Plot a graph using colors based on a partition
function plot_partitioned_graph(nBuses, edgelist, partition, directed=false)
	# Create new colored graph for plotting
	cg = inclist(ExVertex, ExEdge{ExVertex}, is_directed=directed);
	verts = ExVertex[];
	for vid in 1:nBuses
		v = ExVertex(vid, string(vid));
		v.attributes["color"] = COLORS[partition[vid]];
		add_vertex!(cg, v);
		push!(verts, v);
	end
	for (i, edge) in enumerate(edgelist)
		e = ExEdge(i, verts[edge[1]], verts[edge[2]]);
		if partition[edge[1]] == partition[edge[2]]
			e.attributes["color"] = COLORS[partition[edge[1]]];
		else
			e.attributes["color"] = "black";
		end		
		add_edge!(cg, e);
	end
	plot_args(cg);
end

# Plot a graph
function plot_graph(nBuses, edgelist, directed=false)
	# Create new colored graph for plotting
	cg = inclist(ExVertex, ExEdge{ExVertex}, is_directed=directed);
	verts = ExVertex[];
	for vid in 1:nBuses
		v = ExVertex(vid, string(vid));
		add_vertex!(cg, v);
		push!(verts, v);
	end
	for (i, edge) in enumerate(edgelist)
		e = ExEdge(i, verts[edge[1]], verts[edge[2]]);
		add_edge!(cg, e);
	end
	plot_args(cg);
end

# Creates a graph of the buses and branches
function graph_from_data(nBuses, flowin_bus_id, flowout_bus_id, directed=false)
	vertices = 1:nBuses;
	edgelist = collect(zip(flowin_bus_id, flowout_bus_id));

	g = simple_adjlist(nBuses, is_directed=directed);
	for edge in edgelist
		add_edge!(g, edge[1], edge[2]);
	end

	return g, edgelist;
end

# Cluster a graph
function cluster(graph, weights, k, cluster_type)
	if cluster_type == :spectral
		println("TODO: spectral");
	elseif cluster_type == :metis_kway
		objval, part = partGraphKway(graph, k, weights);
		return part;
	elseif cluster_type == :metis_recursive
		objval, part = partGraphRecursive(graph, k, weights);
		return part;
	else
		println("Invalid cluster type $cluster_type");
	end
end
