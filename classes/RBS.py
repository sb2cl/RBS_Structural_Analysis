import json
import math


class RBS:


    def __init__(self):

        # Load JSON file into a dictionary
        with open("data/emopec_data.json", "r") as json_file:
            self.SEQS = json.load(json_file)

    def _get_expression(self, sd_seq):
        """Get the expression level.

        :param str sd_seq:         SD sequence to get expression for.
        """
        sd_seq = sd_seq.upper().replace('U', 'T')

        try:
            raw_expr = self.SEQS[sd_seq]
        except KeyError:
            raw_expr = 0

        return raw_expr

    def predict_spacing(self, leader, max_spacing=12):
        """
        Return elements of the RBS sequence.

        :param str leader:         Leader sequence.
        :param int max_spacing:    Max spacing to look for.

        :return: upstream, sd, spacing, expr
        """
        _leader = leader.upper().replace('U', 'T')
        expr = 0
        spcg = 0

        for i in range(max_spacing):
            # Get the core RBS sequence, start 2 bases before the CDS
            sd = _leader[-7 - i:-1 - i]
            # Get the expression
            _expr = self._get_expression(sd)
            # If the expression is higher than the previous one, update the max value
            if _expr > expr:
                expr = _expr
                spcg = i

        # Get the fractions
        upstream = leader[:-7 - spcg]
        sd = leader[-7 - spcg:-1 - spcg]
        spacing = leader[-1 - spcg:]

        return upstream, sd, spacing, round(expr, 4)
