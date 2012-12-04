from igraph import Graph
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
    def __init__(self, graphs):
        self.graphs = self.convert_graphs(graphs)
        self.n_graphs = len(graphs)
        self.all_calculations = []


    def convert_graphs(self, graphs):
        ans = []
        from scipy.sparse import csc
        from numpy import asarray
        if type(graphs[0]) == csc.csc_matrix:
            for g in range(len(graphs)):
                ans.append(
                        Graph.Weighted_Adjacency(graphs[g].toarray().tolist()))
        elif type(graphs[0]) == list:
            for g in range(len(graphs)):
                ans.append(Graph.Weighted_Adjacency(graphs[g]))
        return asarray(ans)

    def calculate(self, calculations=None):
        if calculations==None or calculations=='all':
            calculations = []
            for i in dir(self):
                if i.startswith('calculate_'):
                    calculations.append(i[10:])

        for var in calculations:
            print("Calculating "+var)
            if type(var)==tuple:
                setattr(self, var[0], var[1](self.graphs))
            elif var.startswith('calculate_'):
                setattr(self, var[10:], getattr(self, var)())
            else:
                setattr(self, var, getattr(self, "calculate_"+var)())

            if var not in self.all_calculations:
                self.all_calculations.append(var)

    def save_to_mat(self, filename):
        from scipy.io import savemat
        data = {}
        for var in self.all_calculations:
            data[var] = getattr(self, var)

        savemat(filename, data)

    def calculator(fn):
        def magic(self):
            from numpy import asarray
            ans = []
            for i in self.graphs:
                ans.append(fn(i))
            return asarray(ans)
        return magic

    @calculator
    def calculate_mincut_value(g):
        return g.mincut_value(capacity=g.es["weight"])

    @calculator
    def calculate_n_infomap(g):
        infomap = g.community_infomap(edge_weights=g.es["weight"])
        return infomap.cluster_graph().vcount()

    @calculator
    def calculate_q_infomap(g):
        infomap = g.community_infomap(edge_weights=g.es["weight"])
        return infomap.q

    @calculator
    def calculate_codelength(g):
        infomap = g.community_infomap(edge_weights=g.es["weight"])
        return infomap.codelength

    @calculator
    def calculate_mean_path_length(g):
        return mean(g.shortest_paths(weights='weight'))

    @calculator
    def calculate_mean_clustering(g):
        return g.transitivity_avglocal_undirected(weights='weight')

    @calculator
    def calculate_rc_itl_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_total_local')
    @calculator
    def calculate_rc_itl_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_total_local')
    @calculator
    def calculate_rc_iNl_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_topN_local')
    @calculator
    def calculate_rc_iNl_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_topN_local')
    @calculator
    def calculate_rc_iNpl_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_topNp_local')
    @calculator
    def calculate_rc_iNpl_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_topNp_local')
    @calculator
    def calculate_rc_iNpwml_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_topNpweightmax_local')
    @calculator
    def calculate_rc_iNpwml_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_topNpweightmax_local')

    @calculator
    def calculate_rc_itg_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_total_global')
    @calculator
    def calculate_rc_itg_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_total_global')
    @calculator
    def calculate_rc_iNg_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_topN_global')
    @calculator
    def calculate_rc_iNg_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_topN_global')
    @calculator
    def calculate_rc_iNpg_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_topNp_global')
    @calculator
    def calculate_rc_iNpg_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_topNp_global')
    @calculator
    def calculate_rc_iNpwmg_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='intensity_topNpweightmax_global')
    @calculator
    def calculate_rc_iNpwmg_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='intensity_topNpweightmax_global')

    @calculator
    def calculate_rc_clustering_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='clustering')
    @calculator
    def calculate_rc_clustering_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='clustering')
    @calculator
    def calculate_rc_n_infomap_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='n_infomap')
    @calculator
    def calculate_rc_n_infomap_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='n_infomap')
    @calculator
    def calculate_rc_q_infomap_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='q_infomap')
    @calculator
    def calculate_rc_q_infomap_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='q_infomap')
    @calculator
    def calculate_rc_codelength_out(g):
        return richclub.rich_club_coefficient(
            g, richness='out_strength',
            club_property='codelength')
    @calculator
    def calculate_rc_codelength_in(g):
        return richclub.rich_club_coefficient(
            g, richness='in_strength',
            club_property='codelength')
    @calculator
    def calculate_rc_in_out_kendalltau(g):
        from scipy.stats import kendalltau
        return  kendalltau(richclub.richness_scores(g, richness='in_strength'),
            richclub.richness_scores(g, richness='out_strength'))

    def calculator_from_first(fn):
        def magic(self):
            from numpy import zeros
            ans = zeros(self.n_graphs)
            first_value = None
            for i in range(1, self.n_graphs):
                if first_value is None:
                    first_value = fn(self.graphs[i])
                else:
                    ans[i] = fn(self.graphs[i], first_value)
            return ans
        return magic

    @calculator_from_first
    def calculate_betweeness_change_from_first(g, first_value=None):
        betweeness_sequence = g.edge_betweenness(weights='weight')
        if first_value:
           from scipy.stats import kendalltau
           return kendalltau(betweeness_sequence, first_value)[0]
        else:
           return betweeness_sequence

    def calculator_derivative(fn):
        def magic(self):
            from numpy import zeros
            ans = zeros(self.n_graphs)
            last_value = None
            for i in range(1, self.n_graphs):
                if last_value is None:
                    last_value = fn(self.graphs[i])
                else:
                    last_value, ans[i] = fn(self.graphs[i], last_value)
            return ans
        return magic

    @calculator_derivative
    def calculate_betweeness_change_kendall(g, last_value=None):
        betweeness_sequence = g.edge_betweenness(weights='weight')
        if last_value:
           from scipy.stats import kendalltau
           return betweeness_sequence, kendalltau(betweeness_sequence, last_value)[0]
        else:
           return betweeness_sequence

    @calculator_derivative
    def calculate_betweeness_change_spearmanr(g, last_value=None):
        betweeness_sequence = g.edge_betweenness(weights='weight')
        if last_value:
           from scipy.stats import spearmanr
           return betweeness_sequence, spearmanr(betweeness_sequence, last_value)[0]
        else:
           return betweeness_sequence

    @calculator_derivative
    def calculate_richness_in_change(g, last_value=None):
        richness_scores = richclub.richness_scores(g, richness='in_strength')
        if last_value:
           from scipy.stats import kendalltau
           return richness_scores, kendalltau(richness_scores, last_value)[0]
        else:
           return richness_scores

    @calculator_derivative
    def calculate_richness_out_change(g, last_value=None):
        richness_scores = richclub.richness_scores(g, richness='out_strength')
        if last_value:
           from scipy.stats import kendalltau
           return richness_scores, kendalltau(richness_scores, last_value)[0]
        else:
           return richness_scores

    @calculator
    def calculate_diameter(g):
        return  g.diameter(weights=g.es["weight"])

    @calculator
    def calculate_diameter_unweighted(g):
        return  g.diameter()

#    def calculate_in_out_strength_cross_correlation(graphs):
#        n_nodes = len(graphs[0].vs)
#        n_generations = len(graphs)
#        from numpy import zeros
#        in_strengths = zeros(n_nodes*n_generations)
#        out_strengths = zeros(n_nodes*n_generations)
#
#        for i in range(n_generations):
#            g = graphs[i]
#            s = g.strength(g.vs, mode=2, weights=g.es["weight"])
#            in_strengths[(i*n_nodes):(n_nodes*(i+1))] =s
#            s = g.strength(g.vs, mode=1, weights=g.es["weight"])
#            out_strengths[(i*n_nodes):(n_nodes*(i+1))] =s
