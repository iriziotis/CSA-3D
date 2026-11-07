Clustering by Geometry
======================

Group catalytic sites by structural similarity.

Hierarchical Clustering
-----------------------

Cluster sites using average-linkage hierarchical clustering::

    sites = list(entry.get_pdbsites(sane_only=True))
    matrix = entry.rmsd_matrix(sites)
    clusters = entry.clustering(matrix)
    
    for cluster_id, site_ids in sorted(clusters.items()):
        print(f"Cluster {cluster_id}: {len(site_ids)} sites")
        for site_id in site_ids:
            print(f"  {site_id}")

Cut Height
----------

Control where to cut the dendrogram::

    # Cut at 70% of tree height (default)
    clusters = entry.clustering(matrix)
    
    # Cut at 50% height (fewer, larger clusters)
    clusters = entry.clustering(matrix, height=0.5)
    
    # Cut at 80% height (more, smaller clusters)
    clusters = entry.clustering(matrix, height=0.8)

Bayesian Clustering
-------------------

Let the algorithm choose the optimal number of clusters::

    clusters = entry.clustering_bayesian(matrix)
    
    # Adjust regularization parameter
    clusters = entry.clustering_bayesian(matrix, l=2.0)

Dendrogram Visualization
------------------------

Plot the dendrogram::

    clusters, linkage_matrix = entry.clustering(
        matrix,
        get_linkage_matrix=True
    )
    
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    
    dendrogram(linkage_matrix, labels=matrix.index)
    plt.show()

Cluster Analysis
----------------

Analyze cluster characteristics::

    for cluster_id, site_ids in clusters.items():
        print(f"Cluster {cluster_id}:")
        
        # Cluster size
        print(f"  Size: {len(site_ids)}")
        
        # Get sites in cluster
        cluster_sites = [entry.get_pdbsite(sid) for sid in site_ids]
        
        # Intra-cluster RMSD
        cluster_matrix = matrix.loc[site_ids, site_ids]
        print(f"  Mean RMSD: {cluster_matrix.values.mean():.2f} Ų")
        
        # Organisms
        organisms = set(s.organism_name for s in cluster_sites)
        print(f"  Organisms: {', '.join(organisms)}")

Working with Clusters
---------------------

Build templates or compare clusters::

    # Get sites from specific cluster
    cluster_sites = [entry.get_pdbsite(sid) for sid in clusters[0]]
    
    # Build template from cluster
    template = entry.create_template(subset=clusters[0])
    
    # Compare clusters structurally
    ref_site = entry.get_pdbsite(clusters[0][0])
    for site_id in clusters[1]:
        site = entry.get_pdbsite(site_id)
        rot, tran, rms, rms_all = ref_site.fit(site)
        print(f"{site_id}: {rms:.2f} Ų inter-cluster")

See Also
--------

- :meth:`csa3d.Entry.clustering` for clustering API
- :meth:`csa3d.Entry.clustering_bayesian` for automatic cluster selection
- :doc:`templates` for building consensus structures
