import json
import math


class RBS:


    def __init__(self):

        # Load JSON file into a dictionary
        with open("data/emopec_data.json", "r") as json_file:
            self.SEQS = json.load(json_file)

    def _spacing_penalty(self, rbs_calc_spacing, optimal_spacing=5,
                        push=(12.2, 2.5, 2.0, 3.0),
                        pull=(0.048, 0.24, 0.0)):
        """Calculate spacing penalty using RBS Calculator formula.

        :param int optimal_spacing: Number of nucleotides which is the opt. spacing.
        :param tuple push: list of c1, c2, c3, c4 if spacing is less than opt.
        :param tuple pull: list of c1, c2, c3 if spacing in more than opt.

        Note that spacing is in RBS calculator spacing, which is one base smaller
        than EMOPEC spacing.
        """
        if rbs_calc_spacing == optimal_spacing:
            return 0.
        elif rbs_calc_spacing < optimal_spacing:
            ds = rbs_calc_spacing - optimal_spacing
            c1, c2, c3, c4 = push
            return c1 / (1.0 + math.exp(c2 * (ds + c3))) ** c4
        else:
            ds = rbs_calc_spacing - optimal_spacing
            c1, c2, c3 = pull
            return c1 * ds ** 2 + c2 * ds + c3

    def _get_expression(self, sd_seq, sd_dist=6, penalty=0.235):
        """Get the expression level including penalty.

        :param str sd_seq:         SD sequence to get expression for.
        :param int sd_dist:        Length of spacing between SD and CDS.
        :param float penalty:      Penalty constant.
        """
        sd_seq = sd_seq.upper().replace('U', 'T')

        try:
            raw_expr = self.SEQS[sd_seq]
        except KeyError:
            raw_expr = 0

        # rbs_calc_spacing
        rcs = sd_dist - 1

        return 10 ** (raw_expr - self._spacing_penalty(rcs) * penalty)

    def predict_spacing(self, leader, penalty=0.235, max_spacing=12):
        """Predict spacing based on a spacing penalty and EMOPEC signal.
        """
        _leader = leader.upper().replace('U', 'T')
        expr = 0
        spcg = 0

        for i in range(max_spacing):
            # Get the core RBS sequence, start 2 bases before the CDS
            sd = _leader[-7 - i:-1 - i]
            # Get the expression
            _expr = self._get_expression(sd, i + 1, penalty=penalty)
            # If the expression is higher than the previous one, update the max value
            if _expr > expr:
                expr = _expr
                spcg = i

        # Get the fractions
        upstream = leader[:-7 - spcg]
        sd = leader[-7 - spcg:-1 - spcg]
        spacing = leader[-1 - spcg:]

        return upstream, sd, spacing, round(expr, 4)
