FILES = {
    '20181015_GSM3258348_control_CD19_B': 'B',
    'GSM3258346_CD_19': 'B',
    'GSM3258346_CD19_B': 'B',
    'GSM3258345_HLA_DR': 'M+DC',
    'GSM3258347_control_HLA_DR': 'M+DC',
    '20170915_GSM2773408_Monocytes_d1': 'MC',
    '20170915_GSM2773409_Monocytes_d2': 'MC',
    '20190108_NK': 'NK',
    '20190108_iNKT': 'T',
    '20190108_MAIT': 'T',
    '20190131_GSM3478792_P5_T': 'T',
    '20181107_GSM3430548_Donor1_ IL-10-producing_Foxp3-_CD4+_T': 'T',
    '20190108_CD4_T': 'T',
    '20190620_GSM3209407_CD4+T': 'T',
    '20190620_GSM3209408_CD4+T_C_CCR5+CD69-': 'T',
    '20190108_CD8': 'T',
    '20180725_GSM3087629_CD8+T_methanol_SSC': 'T',
    '20190108_Vd1': 'T',
    '20190108_Vd2': 'T',
    'file': 'type',
}

with open('../data/interim/geni_shuffled.csv', 'r') as fajlhendel, \
        open('../data/interim/novi.csv', 'w') as autput:
    lines = fajlhendel.readlines()
    for line in lines:
        items = line.split(',')
        try:
            cell_type = FILES['_'.join(items[1].split('_')[1:-2])]
        except:
            cell_type = FILES[items[1]]
        items.insert(2, cell_type)
        autput.write(','.join(items))
