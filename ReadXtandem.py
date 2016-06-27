import os
from pyteomics import tandem
import xml.etree.ElementTree as ET

input_dir = "input"
output_dir = "output"
param_file = "temp.xml"
temaplate_param_file = "xtandem_pref.xml"
coverage_threshold = 0.15

def main():

    for file in os.listdir(output_dir):
        path = output_dir + "/" + file
        os.remove(path)

    for file in os.listdir(input_dir):
        input = input_dir + "/" + file
        output = output_dir + "/" + file
        writeXTandemPrefrences(input, output)
        os.system("./tandem " + param_file)

    for file in os.listdir(output_dir):
        path = output_dir + "/" + file
        xtandemResult = XTandemOutput(path)
        xtandemResult.print()
        print()

def writeXTandemPrefrences(input, output):
    tree = ET.parse(temaplate_param_file)
    root = tree.getroot()
    for param in root:
        if param.attrib["label"] == "output, path":
            param.text = output
        elif param.attrib["label"] == "spectrum, path":
            param.text = input
    tree.write(param_file)

class XTandemOutput:
    path = None
    proteins = {}

    def __init__(self, path):
        self.path = path
        for unit in tandem.read(path):
            for protein in unit["protein"]:
                label = protein["label"]
                peptide = protein["peptide"]
                start = peptide["start"]
                end = peptide["end"]
                if label in self.proteins.keys():
                    self.proteins[label].add(PeptideHit(start, end))
                else:
                    hit = ProteinHit(peptide["peptide"], label)
                    hit.add(PeptideHit(start, end))
                    self.proteins[label] = hit

    def print(self):
        print(self.path)
        i = 0
        for protein in sorted(self.proteins.values(), key=lambda pro: len(pro.hits), reverse=True):
            if (protein.getCoverage() > coverage_threshold):
                protein.print()
            i += 1
            if i == 10: break

class PeptideHit:
    start = 0
    end = 0

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def length(self):
        return self.end - self.start

class ProteinHit:
    hits = []
    protein = None
    label = None
    coverage = []

    def __init__(self, protein, label):
        self.setProtein(protein)
        self.label = label

    def add(self, hit):
        self.hits.append(hit)
        for i in range(hit.start, hit.end):
            self.coverage[i] = True

    def setProtein(self, peptide):
        self.peptide = ""
        self.coverage = []
        for c in peptide:
            if c.isalpha():
                self.peptide += c
                self.coverage.append(False)

    def getCoverage(self):
        coverageCount = 0
        for isCovered in self.coverage:
            if isCovered: coverageCount += 1
        return coverageCount / len(self.coverage)

    def getPeptide(self, hit):
        return self.peptide[hit.start - 1 : hit.end]

    def print(self):
        print("%5.1f%% %-5i %s" % (self.getCoverage() * 100, len(self.hits), self.label))

main()
