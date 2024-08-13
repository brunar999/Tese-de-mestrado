#Python- remove otus with less than 3 reads in 3 replicas of one sample
import pandas as pd
import csv

readfile = "otutab_sing_editado.csv"
outfile = "output_sing_edited_mais_que_2_zero.csv"

wfh = open(outfile, "w")

with open(readfile, "r") as fh:
    reader = csv.DictReader(fh)
    i = 0
    # REORGANIZAR COLUNAS AQUI
    fieldnames = ["ï»¿#OTU ID", "SampleName10-PAEinf", "SampleName16-PAEinf", "SampleName24-PAEinf",
                  "SampleName38-PAEinf",
                  "SampleName11-TGctrl",
                  "SampleName15-TGctrl", "SampleName2-TGctrl", "SampleName3-TGctrl", "SampleName12-Pginf",
                  "SampleName19-Pginf",
                  "SampleName1-Pginf", "SampleName5-Pginf", "SampleName13-PGDctrl", "SampleName17-PGDctrl",
                  "SampleName23-PGDctrl",
                  "SampleName4-PGDctrl", "SampleName14-PGctrl", "SampleName28-PGctrl", "SampleName29-PGctrl",
                  "SampleName39-PGctrl",
                  "SampleName18-Tginf", "SampleName20-Tginf", "SampleName7-Tginf", "SampleName8-Tginf",
                  "SampleName21-Coolctrl",
                  "SampleName26-Coolctrl", "SampleName31-Coolctrl", "SampleName33-Coolctrl",
                  "SampleName22-PGDinf",
                  "SampleName37-PGDinf", "SampleName40-PGDinf", "SampleName6-PGDinf", "SampleName34-PAEctrl",
                  "SampleName35-PAEctrl",
                  "SampleName36-PAEctrl", "SampleName9-PAEctrl", "SampleName25-Coolinf", "SampleName27-Coolinf",
                  "SampleName30-Coolinf", "SampleName32-Coolinf"]
    writer = csv.DictWriter(wfh, fieldnames=fieldnames)
    writer.writeheader()
    for row in reader:
        PAEinf = 0
        TGctrl = 0
        Pginf = 0
        PGDctrl = 0
        PGctrl = 0
        Tginf = 0
        Coolctrl = 0
        PGDinf = 0
        PAEctrl = 0
        Coolinf = 0
        Sample0 = row["ï»¿#OTU ID"]
        Sample1 = row["SampleName1-Pginf"]
        Sample2 = row["SampleName2-TGctrl"]
        Sample3 = row["SampleName3-TGctrl"]
        Sample4 = row["SampleName4-PGDctrl"]
        Sample5 = row["SampleName5-Pginf"]
        Sample6 = row["SampleName6-PGDinf"]
        Sample7 = row["SampleName7-Tginf"]
        Sample8 = row["SampleName8-Tginf"]
        Sample9 = row["SampleName9-PAEctrl"]
        Sample10 = row["SampleName10-PAEinf"]
        Sample11 = row["SampleName11-TGctrl"]
        Sample12 = row["SampleName12-Pginf"]
        Sample13 = row["SampleName13-PGDctrl"]
        Sample14 = row["SampleName14-PGctrl"]
        Sample15 = row["SampleName15-TGctrl"]
        Sample16 = row["SampleName16-PAEinf"]
        Sample17 = row["SampleName17-PGDctrl"]
        Sample18 = row["SampleName18-Tginf"]
        Sample19 = row["SampleName19-Pginf"]
        Sample20 = row["SampleName20-Tginf"]
        Sample21 = row["SampleName21-Coolctrl"]
        Sample22 = row["SampleName22-PGDinf"]
        Sample23 = row["SampleName23-PGDctrl"]
        Sample24 = row["SampleName24-PAEinf"]
        Sample25 = row["SampleName25-Coolinf"]
        Sample26 = row["SampleName26-Coolctrl"]
        Sample27 = row["SampleName27-Coolinf"]
        Sample28 = row["SampleName28-PGctrl"]
        Sample29 = row["SampleName29-PGctrl"]
        Sample30 = row["SampleName30-Coolinf"]
        Sample31 = row["SampleName31-Coolctrl"]
        Sample32 = row["SampleName32-Coolinf"]
        Sample33 = row["SampleName33-Coolctrl"]
        Sample34 = row["SampleName34-PAEctrl"]
        Sample35 = row["SampleName35-PAEctrl"]
        Sample36 = row["SampleName36-PAEctrl"]
        Sample37 = row["SampleName37-PGDinf"]
        Sample38 = row["SampleName38-PAEinf"]
        Sample39 = row["SampleName39-PGctrl"]
        Sample40 = row["SampleName40-PGDinf"]
        # PAEinf
        if Sample10 == "0":
            PAEinf += 1
        if Sample16 == "0":
            PAEinf += 1
        if Sample24 == "0":
            PAEinf += 1
        if Sample38 == "0":
            PAEinf += 1
        # TGctrl
        if Sample11 == "0":
            TGctrl += 1
        if Sample15 == "0":
            TGctrl += 1
        if Sample2 == "0":
            TGctrl += 1
        if Sample3 == "0":
            TGctrl += 1
        # Pginf
        if Sample12 == "0":
            Pginf += 1
        if Sample19 == "0":
            Pginf += 1
        if Sample1 == "0":
            Pginf += 1
        if Sample5 == "0":
            Pginf += 1
        # PGDctrl
        if Sample13 == "0":
            PGDctrl += 1
        if Sample17 == "0":
            PGDctrl += 1
        if Sample23 == "0":
            PGDctrl += 1
        if Sample4 == "0":
            PGDctrl += 1
        # PGctrl
        if Sample14 == "0":
            PGctrl += 1
        if Sample28 == "0":
            PGctrl += 1
        if Sample29 == "0":
            PGctrl += 1
        if Sample39 == "0":
            PGctrl += 1
        # Tginf
        if Sample18 == "0":
            Tginf += 1
        if Sample20 == "0":
            Tginf += 1
        if Sample7 == "0":
            Tginf += 1
        if Sample8 == "0":
            Tginf += 1
        # Coolctrl
        if Sample21 == "0":
            Coolctrl += 1
        if Sample26 == "0":
            Coolctrl += 1
        if Sample31 == "0":
            Coolctrl += 1
        if Sample33 == "0":
            Coolctrl += 1
        # PGDinf
        if Sample22 == "0":
            PGDinf += 1
        if Sample37 == "0":
            PGDinf += 1
        if Sample40 == "0":
            PGDinf += 1
        if Sample6 == "0":
            PGDinf += 1
        # PAEctrl
        if Sample34 == "0":
            PAEctrl += 1
        if Sample35 == "0":
            PAEctrl += 1
        if Sample36 == "0":
            PAEctrl += 1
        if Sample9 == "0":
            PAEctrl += 1
        # Coolinf
        if Sample25 == "0":
            Coolinf += 1
        if Sample27 == "0":
            Coolinf += 1
        if Sample30 == "0":
            Coolinf += 1
        if Sample32 == "0":
            Coolinf += 1
        # Verificação se tem 2 zeros
        if PAEinf > 2:
            Sample10 = "0"
            Sample16 = "0"
            Sample24 = "0"
            Sample38 = "0"
        if TGctrl > 2:
            Sample11 = "0"
            Sample15 = "0"
            Sample2 = "0"
            Sample3 = "0"
        if Pginf > 2:
            Sample12 = "0"
            Sample19 = "0"
            Sample1 = "0"
            Sample5 = "0"
        if PGDctrl > 2:
            Sample13 = "0"
            Sample17 = "0"
            Sample23 = "0"
            Sample4 = "0"
        if PGctrl > 2:
            Sample14 = "0"
            Sample28 = "0"
            Sample29 = "0"
            Sample39 = "0"
        if Tginf > 2:
            Sample18 = "0"
            Sample20 = "0"
            Sample7 = "0"
            Sample8 = "0"
        if Coolctrl > 2:
            Sample21 = "0"
            Sample26 = "0"
            Sample31 = "0"
            Sample33 = "0"
        if PGDinf > 2:
            Sample22 = "0"
            Sample37 = "0"
            Sample40 = "0"
            Sample6 = "0"
        if PAEctrl > 2:
            Sample34 = "0"
            Sample35 = "0"
            Sample36 = "0"
            Sample9 = "0"
        if Coolinf > 2:
            Sample25 = "0"
            Sample27 = "0"
            Sample30 = "0"
            Sample32 = "0"
        # Escrever a OTU
        # REORGANIZAR COLUNAS AQUI
        writer.writerow({"ï»¿#OTU ID":Sample0, "SampleName10-PAEinf":Sample10, "SampleName16-PAEinf":Sample16, "SampleName24-PAEinf":Sample24, "SampleName38-PAEinf":Sample38,
                          "SampleName11-TGctrl":Sample11,
                          "SampleName15-TGctrl":Sample15, "SampleName2-TGctrl":Sample2, "SampleName3-TGctrl":Sample3, "SampleName12-Pginf":Sample12,
                          "SampleName19-Pginf":Sample19,
                          "SampleName1-Pginf":Sample1, "SampleName5-Pginf":Sample5, "SampleName13-PGDctrl":Sample13, "SampleName17-PGDctrl":Sample17,
                          "SampleName23-PGDctrl":Sample23,
                          "SampleName4-PGDctrl":Sample4, "SampleName14-PGctrl":Sample14, "SampleName28-PGctrl":Sample28, "SampleName29-PGctrl":Sample29,
                          "SampleName39-PGctrl":Sample39,
                          "SampleName18-Tginf":Sample18, "SampleName20-Tginf":Sample20, "SampleName7-Tginf":Sample7, "SampleName8-Tginf":Sample8,
                          "SampleName21-Coolctrl":Sample21,
                          "SampleName26-Coolctrl":Sample26, "SampleName31-Coolctrl":Sample31, "SampleName33-Coolctrl":Sample33,
                          "SampleName22-PGDinf":Sample22,
                          "SampleName37-PGDinf":Sample37, "SampleName40-PGDinf":Sample40, "SampleName6-PGDinf":Sample6, "SampleName34-PAEctrl":Sample34,
                          "SampleName35-PAEctrl":Sample35,
                          "SampleName36-PAEctrl":Sample36, "SampleName9-PAEctrl":Sample9, "SampleName25-Coolinf":Sample25, "SampleName27-Coolinf":Sample27,
                          "SampleName30-Coolinf":Sample30, "SampleName32-Coolinf":Sample32})
wfh.close()
