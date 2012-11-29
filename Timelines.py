from igraph import Graph
from scipy.stats import kendalltau, spearmanr
from numpy import mean, std, empty, shape
import richclub


class Timelines:
    def __init__(self, with_control=True):
        self.tls = []
        self.with_control = with_control

        if self.with_control:
            self.controls = []

    def add_timeline(self, tl, control_tls=None):
        self.tls.append(tl)
        if self.with_control:
            self.controls.append(control_tls)

    def add_timelines(self, tl, control_tls=None):
        for i in tl:
            self.add_timeline(i)
        if self.with_control:
            for i in control_tls:
                self.controls.append(i)

    def raw_data(self, name):
        n_tls = len(self.tls)
        raw_data = empty((n_tls,)+shape(getattr(self.tls[0], name)))

        for i in range(n_tls):
            raw_data[i] = getattr(self.tls[i], name)

        return raw_data

    def raw_control(self, name):
        n_tls = len(self.tls)
        raw_control = empty((n_tls,)+
            shape(self.controls[0].raw_data(name)))

        for i in range(n_tls):
            raw_control[i] = self.controls[i].raw_data(name)

        return raw_control

    def mean_control(self, name):
        raw_control = self.raw_control(name)
        mean_control = mean(raw_control, axis=1)

        return mean_control


    def data(self, name):
        raw_data = self.raw_data(name)
        m = mean(raw_data, axis=0)
        s = std(raw_data, axis=0)

        return m, s

    def control(self, name):
        mean_control = self.mean_control(name)
        m = mean(mean_control, axis=0)
        s = std(mean_control, axis=0)

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

        self.rc_itl_out = []
        self.rc_itl_in = []
        self.rc_iNl_out = []
        self.rc_iNl_in = []
        self.rc_iNpl_out = []
        self.rc_iNpl_in = []
        self.rc_iNpwml_out = []
        self.rc_iNpwml_in = []
        self.rc_itg_out = []
        self.rc_itg_in = []
        self.rc_iNg_out = []
        self.rc_iNg_in = []
        self.rc_iNpg_out = []
        self.rc_iNpg_in = []
        self.rc_iNpwmg_out = []
        self.rc_iNpwmg_in = []
        self.rc_clustering_out = []
        self.rc_clustering_in = []
        self.rc_n_infomap_out = []
        self.rc_n_infomap_in = []
        self.rc_q_infomap_out = []
        self.rc_q_infomap_in = []
        self.rc_codelength_out = []
        self.rc_codelength_in = []

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
        self.mean_clustering.append(
                g.transitivity_avglocal_undirected(weights='weight'))

        betweeness_sequence = g.edge_betweenness(weights='weight')

        if self.last_betweeness is None:
            self.last_betweeness = betweeness_sequence
        else:
            self.betweeness_change_kendall.append(
                    kendalltau(self.last_betweeness, betweeness_sequence)[0])
            self.betweeness_change_spearmanr.append(
                    spearmanr(self.last_betweeness, betweeness_sequence)[0])

        self.rc_itl_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_total_local'))
        self.rc_itl_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_total_local'))
        self.rc_iNl_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_topN_local'))
        self.rc_iNl_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_topN_local'))
        self.rc_iNpl_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_topNp_local'))
        self.rc_iNpl_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_topNp_local'))
        self.rc_iNpwml_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_topNpweightmax_local'))
        self.rc_iNpwml_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_topNpweightmax_local'))

        self.rc_itg_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_total_global'))
        self.rc_itg_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_total_global'))
        self.rc_iNg_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_topN_global'))
        self.rc_iNg_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_topN_global'))
        self.rc_iNpg_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_topNp_global'))
        self.rc_iNpg_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_topNp_global'))
        self.rc_iNpwmg_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='intensity_topNpweightmax_global'))
        self.rc_iNpwmg_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='intensity_topNpweightmax_global'))

        self.rc_clustering_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='clustering'))
        self.rc_clustering_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='clustering'))
        self.rc_n_infomap_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='n_infomap'))
        self.rc_n_infomap_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='n_infomap'))
        self.rc_q_infomap_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='q_infomap'))
        self.rc_q_infomap_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='q_infomap'))
        self.rc_codelength_out.append(
                richclub.rich_club_coefficient(
                    g, richness='out_strength',
                    club_property='codelength'))
        self.rc_codelength_in.append(
                richclub.rich_club_coefficient(
                    g, richness='in_strength',
                    club_property='codelength'))

    def close_gens(self):
        from numpy import asarray
        self.n_infomap = asarray(self.n_infomap)
        self.q_infomap = asarray(self.q_infomap)
        self.codelength = asarray(self.codelength)
        self.mean_path_length = asarray(self.mean_path_length)
        self.mean_clustering = asarray(self.mean_clustering)
        self.betweeness_change_kendall = asarray(self.betweeness_change_kendall)
        self.betweeness_change_spearmanr = asarray(self.betweeness_change_spearmanr)

        self.rc_itl_out = asarray(self.rc_itl_out)
        self.rc_itl_in = asarray(self.rc_itl_in)
        self.rc_iNl_out = asarray(self.rc_iNl_out)
        self.rc_iNl_in = asarray(self.rc_iNl_in)
        self.rc_iNpl_out = asarray(self.rc_iNpl_out)
        self.rc_iNpl_in = asarray(self.rc_iNpl_in)
        self.rc_iNpwml_out = asarray(self.rc_iNpwml_out)
        self.rc_iNpwml_in = asarray(self.rc_iNpwml_in)
        self.rc_itg_out = asarray(self.rc_itg_out)
        self.rc_itg_in = asarray(self.rc_itg_in)
        self.rc_iNg_out = asarray(self.rc_iNg_out)
        self.rc_iNg_in = asarray(self.rc_iNg_in)
        self.rc_iNpg_out = asarray(self.rc_iNpg_out)
        self.rc_iNpg_in = asarray(self.rc_iNpg_in)
        self.rc_iNpwmg_out = asarray(self.rc_iNpwmg_out)
        self.rc_iNpwmg_in = asarray(self.rc_iNpwmg_in)
        self.rc_clustering_out = asarray(self.rc_clustering_out)
        self.rc_clustering_in = asarray(self.rc_clustering_in)
        self.rc_n_infomap_out = asarray(self.rc_n_infomap_out)
        self.rc_n_infomap_in = asarray(self.rc_n_infomap_in)
        self.rc_q_infomap_out = asarray(self.rc_q_infomap_out)
        self.rc_q_infomap_in = asarray(self.rc_q_infomap_in)
        self.rc_codelength_out = asarray(self.rc_codelength_out)
        self.rc_codelength_in = asarray(self.rc_codelength_in)
