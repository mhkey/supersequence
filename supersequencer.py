import itertools


class Node:
    def __init__(self, key: str | int, data=None):
        self.id = key
        self.data = data
        self.start: bool = False
        self.end: bool = False
        self.prev = set()
        self.next = set()

    def __repr__(self):
        return str(self.id)

    def __eq__(self, other):
        return other.id == self.id if isinstance(other, Node) else self.id == other

    def __hash__(self):
        return hash(self.id)


def is_super_sequence(short_seq: list, long_seq: list) -> bool:
    it = iter(long_seq)
    return all(item in it for item in short_seq)


def eliminate_redundant_sequences(sequences: list) -> list:
    sequences.sort(key=len)
    result = []
    for sequence in sequences:
        if not any(is_super_sequence(sequence, other) for other in result):
            result.append(sequence)
    return result


def get_all_nodes(sequences: list) -> list[Node]:
    nodes = {}
    for sequence in sequences:
        for i, node_id in enumerate(sequence):
            node = nodes.setdefault(node_id, Node(node_id))
            if i == 0:
                node.start = True
            if i == len(sequence) - 1:
                node.end = True
            if i > 0:
                prev_node_id = sequence[i - 1]
                node.prev.add(prev_node_id)
                nodes[prev_node_id].next.add(node_id)
    return list(nodes.values())


def get_sections(all_nodes: list[Node], sequences: list) -> list:
    """groups sequences into sub sequence sections based on all nodes start, end nodes flags
    previous, next node counts. """
    sections = []
    for i, seq_section in enumerate(sequences):
        sub_seq = []
        sections.append([])
        for p in seq_section:
            node = all_nodes[all_nodes.index(p)]
            if (len(node.prev) > 1 or node.start or (node.end and node.next)) and sub_seq:
                sections[i].append(sub_seq)
                sub_seq = [node]
            else:
                sub_seq.append(node)
            if (len(node.next) > 1 or node.end) and sub_seq:
                sections[i].append(sub_seq)
                sub_seq = []

    return sorted(sections, key=len, reverse=True)


def shift_orphans(super_seq: list) -> list:
    for _ in range(10):
        shifted = False
        for i in range(1, len(super_seq)):
            if len(super_seq[i]) == 1 and not super_seq[i][0].start:
                prev_node = super_seq[i - 1][-1]
                if prev_node not in super_seq[i][0].prev:
                    super_seq[i - 1], super_seq[i] = super_seq[i], super_seq[i - 1]
                    shifted = True
                    break
        if not shifted:
            break
    return super_seq


def remove_duplicates(sequences: list) -> list:
    to_remove = []  # build list to remove
    for node in sequences:
        multi = [index for index, value in enumerate(sequences) if value == node]
        if len(multi) > 1:
            if all([True
                    if n in sequences[multi[-1]:] and n not in sequences[multi[0]:multi[-1]]
                    else False for n in node.next]) and multi[0] not in to_remove:
                to_remove.append(multi[0])
    # remove
    for index, value in enumerate(to_remove):
        sequences.pop(to_remove[index] - index)
    return sequences


def super_sequence(sequences_list):

    sequences = eliminate_redundant_sequences([[node for node in seq] for seq in sequences_list])

    all_nodes = get_all_nodes(sequences)
    subsequences_list = get_sections(all_nodes, sequences)

    super_sequence = []
    for sequence in subsequences_list:
        i = 0
        tmp_sequences = []
        for sub_sequence in sequence:
            try:
                i = super_sequence.index(sub_sequence, i) + 1
                for sub in tmp_sequences:
                    super_sequence.insert(i - 1, sub)
                    i += 1
                tmp_sequences = []
            except ValueError:
                tmp_sequences.append(sub_sequence)

        for k, v in enumerate(tmp_sequences):
            super_sequence.insert(i + k, v)

    super_sequence = shift_orphans(super_sequence)

    return list(itertools.chain(*super_sequence))
