import itertools
from copy import deepcopy


class Point:
    def __init__(self, seq_point=None):
        if isinstance(seq_point, dict):
            self.id = seq_point.get('tsn', seq_point)  # identifier Name of object
        else:
            self.id = seq_point
        self.data = seq_point  # data of sequence
        self.start = False  # Flag if Point is the First Point of any sub_sequence
        self.end = False  # Flag if Point is the End Point of any sub_sequence
        self.prev = []  # List of Points that preceded current
        self.next = []  # List of Points that follows current

    def __repr__(self):
        return str(self.id)

    def __eq__(self, other):
        return other.id == self.id if isinstance(other, Point) else self.id == other


def supersequence(sequences_list):
    """ Generates common super sequence """

    def is_supersequence(xs, ys):
        """ True if y is supersequence of x """
        idx = 0
        try:
            for item in ys:
                idx = xs.index(item, idx) + 1
        except ValueError:
            return False
        return True

    def eliminate_redundant_sequences(seqs):
        """ Remove Redundant sequences (duplicates and is_super_sequence matches)"""
        redundant_sequences = []
        seqs.sort(key=lambda xs: len(xs))
        for idx, req_red in enumerate(seqs):
            if idx not in redundant_sequences:
                for j in range(idx+1, len(seqs)):
                    if req_red == seqs[j]:
                        redundant_sequences.append(idx)
                        break
                    elif is_supersequence(seqs[j], req_red):
                        redundant_sequences.append(idx)
                        break
        return [vx for ix, vx in enumerate(seqs) if ix not in redundant_sequences]

    def get_allpoints(seqs):
        """ get all points in all sequences flag start\end - list prev and next points in sequence """
        prev_point = None
        points = []
        for seq_ap in seqs:
            for ix, point in enumerate(seq_ap):
                if point not in points:
                    points.append(point)
                    index = len(points)-1
                else:
                    index = points.index(point)

                if ix == 0:
                    points[index].start = True
                    prev_point = None

                if ix == len(seq_ap)-1:
                    points[index].end = True

                if ix > 0:
                    if prev_point is not None:
                        if prev_point not in points[index].prev:
                            points[index].prev.append(prev_point)

                        if point not in points[points.index(prev_point)].next:
                            points[points.index(prev_point)].next.append(point)
                prev_point = point
        return points

    def get_sections(allpoints, seq1):
        """groups sequences into sub sequence sections based on all points start, end points flags
        previous, next point counts. """
        seq_list = []
        for ix, seq_section in enumerate(seq1):
            sub_seq = []
            seq_list.append([])
            for p in seq_section:
                point = allpoints[allpoints.index(p)]
                if (len(point.prev) > 1 or point.start or (point.end and point.next)) and sub_seq:
                    seq_list[ix].append(sub_seq)
                    sub_seq = [point]
                else:
                    sub_seq.append(point)
                if (len(point.next) > 1 or point.end) and sub_seq:
                    seq_list[ix].append(sub_seq)
                    sub_seq = []

        seq_list = sorted(seq_list, key=lambda xx: len(xx), reverse=True)  # Sorts by how many junctions each trip has..
        # return seq_list
        # try:
        #    return sorted(seq_list, key=lambda x: len(x[0][-1].next[0].prev), reverse=False)
        # except:

        return seq_list

    def remove_duplicates(seqs):
        toremove = []  # build list to remove
        for point in seqs:
            multi = [index for index, value in enumerate(seqs) if value == point]
            if len(multi) > 1:
                if all([True
                        if n in seqs[multi[-1]:] and n not in seqs[multi[0]:multi[-1]]
                        else False for n in point.next]) and multi[0] not in toremove:
                    toremove.append(multi[0])
        # remove
        for index, value in enumerate(toremove):
            seqs.pop(toremove[index]-index)
        return seqs

    sequences = eliminate_redundant_sequences([[Point(p) for p in s] for s in sequences_list])

    all_points = get_allpoints(sequences)
    subsequences_list = get_sections(all_points, sequences)

    super_sequence = []
    for sequence in subsequences_list:
        i = 0  # position index
        tmp_sequences = []
        for sub_sequence in sequence:
            try:  # EAFP
                i = super_sequence.index(sub_sequence, i) + 1
                for sub in tmp_sequences:
                    super_sequence.insert(i-1, sub)
                    i += 1
                tmp_sequences = []
            except ValueError:  # [x] not in list
                tmp_sequences.append(sub_sequence)

        # remaining sub sequences insert inline
        for k, v in enumerate(tmp_sequences):
            super_sequence.insert(i+k, v)

    # TODO: Weighting shift single orphaned stops closer to previous sequence segment. What others weightings todo?
    ttt = 10  # number of passes.
    while ttt > 0:
        ttt -= 1
        shift1 = []
        for i, x in enumerate(super_sequence):
            if len(x) == 1 and not x[0].start and i > 0:
                if super_sequence[i-1][-1] not in x[0].prev:
                    shift1.append(i)
                    break
        if shift1:
            for x in shift1:
                super_sequence[x-1], super_sequence[x] = super_sequence[x], super_sequence[x-1]
        else:
            break

    final = remove_duplicates([deepcopy(t) for t in list(itertools.chain(*super_sequence))])
    return final


if __name__ == '__main__':
    seq = [
        ['20006', '20471', '211158', '20462', '20461', '21371', '21122', '21141', '21271', '21151'],
        ['2000442', '21271'],
        ['20006', '20471', '211158', '20462', '20461', '21371', '21122', '21141', '21271'],
        ['20006', '2000442', '20009', '20471', '211158', '20462', '20461', '21371', '21122', '21141', '21271'],
        ['20005', '2000442', '20009', '20471', '211158', '20462', '20461', '21371', '21122', '21141', '21271']

    ]
    print(supersequence(seq))
