from pybiotk.utils import split_discontinuous


wason_crick_pair = {"A": "U", "C": "G", "G": "C", "U": "A",}


class Seed:
    def __init__(self, seed_pos_list, gu_pos_list, gap_pos_list, gap_len_list):
        self.seed = seed_pos_list
        self.gap_pos = gap_pos_list
        self.gap_len = gap_len_list
        self.gu = gu_pos_list
        self.seed_start = seed_pos_list[0]
        self.seed_stop = seed_pos_list[-1]
        self.seed_len = self.seed_stop - self.seed_start + 1
        self.gu_count = len(self.gu)
        self.gap_dict = dict(zip(self.gap_pos, self.gap_len))
        
    def __repr__(self):
        return f"{self.seed_start}-{self.seed_stop}"
        
    __str__ = __repr__
        
    def window_in_seed(self, start, end):
        window_set = set(range(start, end+1))
        seed_set = set(self.seed)
        if window_set - seed_set:
            return False
        else:
            return True
    
    def window_gu_gap(self, start, end):
        window_set = set(range(start, end+1))
        gu_set = set(self.gu)
        gap_set = set(self.gap_pos)
        
        window_gu = list(window_set & gu_set)
        window_gap = list(window_set & gap_set)
        window_gap_len = list(self.gap_dict[i] for i in window_gap)
        
        return window_gu, window_gap, window_gap_len
        
        
class HybridStructure:
    def __init__(self, query_sequence, query_structure, query_position, target_sequence, target_structure):
        self.query_sequence = query_sequence.replace("T", "U")
        self.query_structure = query_structure
        self.query_position = query_position
        self.target_sequence = target_sequence.replace("T", "U")
        self.target_structure = target_structure
        self.query_length = len(self.query_sequence)
        self.target_length = len(self.target_sequence)

    @classmethod
    def get_basepair_pos_list(cls, structure):
        basepairlist = []
        for i, symbol in enumerate(structure):
            if not symbol == '.':
                basepairlist.append(i)
        return basepairlist

    def query_pair_pos(self):
        return self.get_basepair_pos_list(self.query_structure)

    def target_pair_pos(self):
        return list(reversed(self.get_basepair_pos_list(self.target_structure)))

    def classfy_basepair(self):
        wason_crick_pos = []
        gu_pos = []
        query_pos = self.query_pair_pos()
        target_pos = self.target_pair_pos()
        for q_pos, t_pos in zip(query_pos, target_pos):
            q_base = self.query_sequence[q_pos]
            t_base = self.target_sequence[t_pos]
            if wason_crick_pair[q_base] == t_base:
                wason_crick_pos.append((q_pos, t_pos))
            elif (q_base, t_base) == ("G", "U") or (q_base, t_base) == ("U", "G"):
                gu_pos.append((q_pos, t_pos))
            else:
                raise Exception(print(self.query_sequence, self.target_sequence, q_pos, t_pos, q_base, t_base))
        return wason_crick_pos, gu_pos

    def get_wason_crick_pos(self):
        wason_crick_pos, _ = self.classfy_basepair()
        return wason_crick_pos

    def get_gu_pos(self):
        _, gu_pos = self.classfy_basepair()
        return gu_pos
    
    def get_basepair_num(self):
        wason_crick_pos, gu_pos = self.classfy_basepair()
        return (len(wason_crick_pos), len(gu_pos))

    def get_wason_crick_num(self):
        num, _ = self.get_basepair_num()
        return num
    
    def get_gu_num(self):
        _, num = self.get_basepair_num()
        return num

    def structure2vector(self, wason_crick=1, gu=0, notpair=0):
        wason_crick_pos, gu_pos = self.classfy_basepair()
        wason_crick_dict = dict(wason_crick_pos)
        gu_dict = dict(gu_pos)
        query_vector = []
        target_vector = []
        for i in range(len(self.query_sequence)):
            if i in wason_crick_dict:
                query_vector.append(wason_crick)
            elif i in gu_dict:
                query_vector.append(gu)
            else:
                query_vector.append(notpair)

        for i in range(len(self.target_sequence)):
            if i in wason_crick_dict.values():
                target_vector.append(wason_crick)
            elif i in gu_dict.values():
                target_vector.append(gu)
            else:
                target_vector.append(notpair)

        return query_vector, target_vector

    def query_structure2vector(self, wason_crick=1, gu=0.5, notpair=0):
        query_vector, _ = self.structure2vector(
            wason_crick=wason_crick, gu=gu, notpair=notpair)
        return query_vector

    def target_structure2vector(self, wason_crick=1, gu=0.5, notpair=0):
        _, target_vector = self.structure2vector(
            wason_crick=wason_crick, gu=gu, notpair=notpair)
        return target_vector

    def mark_gu_structure(self, query_gu_symbol="*", target_gu_symbol="*"):
        query_structure_list = list(self.query_structure)
        target_structure_list = list(self.target_structure)
        gu_pos = self.get_gu_pos()
        for (q_pos, t_pos) in gu_pos:
            query_structure_list[q_pos] = query_gu_symbol
            target_structure_list[t_pos] = target_gu_symbol
        return "".join(query_structure_list), "".join(target_structure_list)

    def mark_query_gu_structure(self, gu_symbol="*"):
        query_structure, _ = self.mark_gu_structure(query_gu_symbol=gu_symbol)
        return query_structure

    def mark_target_gu_structure(self, gu_symbol="*"):
        _, target_structure = self.mark_gu_structure(target_gu_symbol=gu_symbol)
        return target_structure
    
    def get_structure(self):
        query_structure_gu, target_structure_gu = self.mark_gu_structure(query_gu_symbol="*", target_gu_symbol="*")
        if self.query_position == "left":
            return (f">query&target\n{self.query_sequence}&{self.target_sequence}\n"
                    f"{self.query_structure}&{self.target_structure}\n"
                    f"{query_structure_gu}&{target_structure_gu}\n")
        elif self.query_position == "right":
            return (f">target&query\n{self.target_sequence}&{self.query_sequence}\n"
                    f"{self.target_structure}&{self.query_structure}\n"
                    f"{target_structure_gu}&{query_structure_gu}\n")
        else:
            raise Exception("Wrong query position.")
    
    def query_seed(self, min_szie=4):
        query_pos = self.query_pair_pos()
        target_pos = self.target_pair_pos()
        gu_pos = self.get_gu_pos()
        basepair_dict = dict(zip(query_pos, target_pos))
        gu_dict = dict(gu_pos)
        seed_list = []
        for seed_region in split_discontinuous(query_pos):
            if len(seed_region) < min_szie:
                continue
            gu_seed = []
            gap_seed = []
            gap_length = []
            last_index = []
            for index in seed_region:
                if index in gu_dict:
                    gu_seed.append(index)
                if not last_index:
                    last_index.append(index)
                    continue
                last = last_index[-1]
                distance = abs(basepair_dict[index] - basepair_dict[last]) - 1
                if distance > 0:
                    gap_seed.append(last)
                    gap_length.append(distance)
                last_index.append(index)
            seed = Seed(seed_pos_list=seed_region, gu_pos_list=gu_seed, gap_pos_list=gap_seed, gap_len_list=gap_length)
            seed_list.append(seed)
            
        return seed_list
        
    def query_seed_noGU(self, min_szie=4):
        wason_crick_pos = self.get_wason_crick_pos()
        wason_crick_dict = dict(wason_crick_pos)
        seed_list = []
        for seed_region in split_discontinuous(wason_crick_dict.keys()):
            if len(seed_region) < min_szie:
                continue
            gu_seed = []
            gap_seed = []
            gap_length = []
            last_index = []
            for index in seed_region:
                if not last_index:
                    last_index.append(index)
                    continue
                last = last_index[-1]
                distance = abs(wason_crick_dict[index] - wason_crick_dict[last]) - 1
                if distance > 0:
                    gap_seed.append(last)
                    gap_length.append(distance)
                last_index.append(index)
            seed = Seed(seed_pos_list=seed_region, gu_pos_list=gu_seed, gap_pos_list=gap_seed, gap_len_list=gap_length)
            seed_list.append(seed)
            
        return seed_list
    
    def query_window(self, size=4, max_gu=None, max_gap=None, max_gap_len=None):
        seed_list = self.query_seed()
        window_list = []
        if not seed_list:
            return window_list
        for i in range(self.query_length - size + 1):
            gu = []
            gap = []
            gap_length = []
            for seed in seed_list:
                region_flag = True
                if seed.window_in_seed(i, i+size-1):
                    gu, gap, gap_length = seed.window_gu_gap(i, i+size-1)
                    break
                else:
                    region_flag = False
            if not region_flag:
                continue
            if max_gu is not None and len(gu) > max_gu:
                continue
            if max_gap is not None and len(gap) > max_gap:
                continue
            longest_gap = max(gap_length) if gap_length else 0
            if max_gap_len is not None and longest_gap > max_gap_len:
                continue 
            window_list.append((i, i+size-1))
        return window_list

