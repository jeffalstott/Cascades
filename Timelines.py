from igraph import Graph
from scipy.stats import kendalltau, spearmanr
from numpy import mean, std, asarray
import richclub


class Timelines:
    def __init__(self, with_control=True):
        self.tls = []
        self.control = with_control

        if self.control:
            self.controls = []

    def add_timeline(self, tl, control_tls=None):
        self.tls.append(tl)
        if self.control:
            self.controls.append(control_tls)

    def add_timelines(self, tl, control_tls=None):
        for i in tl:
            self.add_timeline(i)
        if self.control:
            for i in control_tls:
                self.controls.append(i)

    def __getattr__(self, name, version=None):

        if version is None:
            if self.control:
                version = 'normalized'
            else:
                version = 'raw'

        if version in ['normalized', 'raw']:
            raw = []
            for i in self.tls:
                raw.append(getattr(i, name))
        if version in ['normalized', 'control']:
            control = []
            for i in self.controls:
                control.append(getattr(i, name)[0])

        if version == 'raw':
            m = mean(asarray(raw), axis=0)
            s = std(asarray(raw), axis=0)
        elif version == 'control':
            m = mean(asarray(control), axis=0)
            s = std(asarray(control), axis=0)
        elif version == 'normalized':
            data = asarray(raw) / asarray(control)
            m = mean(data, axis=0)
            s = std(data, axis=0)

        return m, s


class Timeline:
    def __init__(self):
        self.q_betweeness = []
        self.n_betweeness = []
        self.q_walktrap = []
        self.n_walktrap = []
        self.q_infomap = []
        self.n_infomap = []
        self.codelength = []

        self.mean_path_length = []
        self.mean_clustering = []

        self.betweeness_change_kendall = []
        self.betweeness_change_spearmanr = []

        self.rc_out = []  # zeros([n_nets, n_runs, 9])
        self.rc_in = []  # zeros([n_nets, n_runs, 9])
        self.rc_out_int = []
        self.rc_in_int = []

        self.last_betweeness = None

    def add_gen(self, g):
        if type(g) == list:
            g = Graph.Weighted_Adjacency(g)
        #b = g.community_edge_betweenness(directed=True, weights=g.es["weight"])
        #n_betweeness.append(b.optimal_count)
        #q_betweeness.append(b.as_clustering().q)

        #w = g.community_walktrap(weights=g.es["weight"], steps=100)
        #n_walktrap.append(w.optimal_count)
        #q_walktrap.append(w.as_clustering().q)

        infomap = g.community_infomap(edge_weights=g.es["weight"])
        self.n_infomap.append(infomap.cluster_graph().vcount())
        self.q_infomap.append(infomap.q)
        self.codelength.append(infomap.codelength)

        self.mean_path_length.append(mean(g.shortest_paths(weights='weight')))
        self.mean_clustering.append(mean(g.transitivity_local_undirected(weights='weight')))

        betweeness_sequence = g.edge_betweenness(weights='weight')

        if self.last_betweeness is None:
            self.last_betweeness = betweeness_sequence
        else:

            #print shape(self.last_betweeness)
            #print shape(betweeness_sequence)

            #kendalltau(self.last_betweeness, betweeness_sequence)

            self.betweeness_change_kendall.append(kendalltau(self.last_betweeness, betweeness_sequence)[0])
            self.betweeness_change_spearmanr.append(spearmanr(self.last_betweeness, betweeness_sequence)[0])

        self.rc_out.append(richclub.rich_club_coefficient(g, scores_name='out_strength', rewire=False))
        self.rc_in.append(richclub.rich_club_coefficient(g, scores_name='in_strength', rewire=False))
        self.rc_out_int.append(sum(self.rc_out[-1] - 1))
        self.rc_in_int.append(sum(self.rc_in[-1] - 1))
