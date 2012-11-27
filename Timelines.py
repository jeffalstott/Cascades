from igraph import Graph
from scipy.stats import kendalltau, spearmanr
from numpy import mean, std, empty
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

    def raw_data(self, name):
        n_tls = len(self.tls)
        raw_data = empty(n_tls,
            dtype=type(getattr(self.tls[0], name)))

        for i in range(n_tls):
            raw_data[i] = getattr(self.tls[i], name)

        return raw_data

    def raw_control(self, name):
        n_tls = len(self.tls)
        raw_control = empty(n_tls,
            dtype=type(self.controls[0].raw(name)))

        for i in range(n_tls):
            raw_control[i] = self.controls[i].raw(name)

        return raw_control

    def data(self, name):
        raw_data = self.raw_data(name)
        m = mean(raw_data, axis=0)
        s = std(raw_data, axis=0)

        return m, s

    def control(self, name):
        raw_control = self.raw_control(name)
        raw_control = mean(raw_control, axis=1)
        m = mean(raw_control, axis=0)
        s = std(raw_control, axis=0)

        return m, s

    def raw_normalized(self, name):

        #each data point should be paired to its own control
        raw_data = self.raw_data(name)
        raw_control = self.raw_control(name)
        raw_control = mean(raw_control, axis=1)
        raw_normalized = raw_data / raw_control

        return raw_normalized

    def normalized(self, name):
        raw_normalized = self.raw_normalized(name)

        m = mean(raw_normalized, axis=0)
        s = std(raw_normalized, axis=0)

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
            self.betweeness_change_kendall.append(kendalltau(self.last_betweeness, betweeness_sequence)[0])
            self.betweeness_change_spearmanr.append(spearmanr(self.last_betweeness, betweeness_sequence)[0])

        self.rc_out.append(richclub.rich_club_coefficient(g, scores_name='out_strength', rewire=False))
        self.rc_in.append(richclub.rich_club_coefficient(g, scores_name='in_strength', rewire=False))
        self.rc_out_int.append(sum(self.rc_out[-1] - 1))
        self.rc_in_int.append(sum(self.rc_in[-1] - 1))

    def close_gens(self):
        from numpy import asarray
        self.n_infomap = asarray(self.n_infomap)
        self.q_infomap = asarray(self.q_infomap)
        self.codelength = asarray(self.codelength)
        self.mean_path_length = asarray(self.mean_path_length)
        self.mean_clustering = asarray(self.mean_clustering)
        self.betweeness_change_kendall = asarray(self.betweeness_change_kendall)
        self.betweeness_change_spearmanr = asarray(self.betweeness_change_spearmanr)
        self.rc_out = asarray(self.rc_out)
        self.rc_in = asarray(self.rc_in)
        self.rc_out_int = asarray(self.rc_out_int)
        self.rc_in_int = asarray(self.rc_in_int)
