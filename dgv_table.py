import os
import database

class dgv_map:
  def __init__(self, dgv):
      dgv = [None if d == 'None' else d for d in dgv]

      self.chr = dgv[0]
      self.start = dgv[1]
      self.end = dgv[2]
      self.state = dgv[3]
      self.id = dgv[4]
      self.type = dgv[5]
      self.num_variants = dgv[6]
      self.num_samples = dgv[7]
      self.num_samples_multicounted = dgv[8]
      self.num_studies = dgv[9]
      self.variants = dgv[10]
      self.samples = dgv[11]
      self.studies = dgv[12]
      self.African = dgv[13]
      self.Asia = dgv[14]
      self.European = dgv[15]
      self.Mexican = dgv[16]
      self.Middle_East = dgv[17]
      self.Native_American = dgv[18]
      self.Oceania = dgv[19]
      self.South_American = dgv[20]

      def __str__(self):
          return ",".join([self.chr, self.start, self.end, self.state, self.id, self.type, self.num_variants,
                self.num_samples,self.num_samples_multicounted, self.num_studies, self.variants, self.samples,
                self.studies, self.African, self.Asia, self.European, self.Mexican, self.Middle_East, self.Native_American,
                self.Oceania, self.South_American])

class cnv_custom_map:
    def __init__(self, cnv_map):
        cnv_map = [None if c == 'None' else c for c in cnv_map]

        self.chr = cnv_map[0]
        self.start = cnv_map[1]
        self.end = cnv_map[2]
        self.opt_field = cnv_map[3]

    def __str__(self):
        return ",".join([self.chr, sefl.start, self.end, self.opt_fieldspi])
