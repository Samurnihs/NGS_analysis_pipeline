import numpy as np
import pandas as pd
import argparse


def isbrit(dfline): # checking of affiliation to UK strain
    if (dfline['del HV69-70'] == 'ДА') and (dfline['del Y144'] == 'ДА') and (dfline['N501Y'] == 'ДА') and (dfline['A570D'] == 'ДА'):
                    return 'ДА'

    elif (dfline['del HV69-70'] == 'НЕТ') or (dfline['del Y144'] == 'НЕТ') or (dfline['N501Y'] == 'НЕТ') or (dfline['A570D'] == 'НЕТ'):
    	return 'НЕТ'
    
    else:
        return 'не определено'


def issar(dfline): # checking of affiliation to SAR strain
    #print(dfline[4], dfline[6], dfline[5])
    if (dfline['D80A'] == 'ДА') and (dfline['N501Y'] == 'ДА') and (dfline['E484K'] == 'ДА'):
                return 'ДА'
    
    elif (dfline['D80A'] == 'НЕТ') or (dfline['N501Y'] == 'НЕТ') or (dfline['E484K'] == 'НЕТ'):
        return 'НЕТ'
    
    else:
        return 'не определено'



def isbras(dfline): # checking of affiliation to Brazil strain
    if (dfline['D138Y'] == 'ДА') and (dfline['N501Y'] == 'ДА') and (dfline['E484K'] == 'ДА'):
    	return 'ДА'
    
    elif (dfline['D138Y'] == 'НЕТ') or (dfline['N501Y'] == 'НЕТ') or (dfline['E484K'] == 'НЕТ'):
        return 'НЕТ'
    
    else:
        return 'не определено'


def isind(dfline):
    if (dfline['L452R'] == 'ДА') and (dfline['T478K'] == 'ДА') and (dfline['P681R'] == 'ДА'):
    	return 'ДА'
    elif (dfline['L452R'] == 'НЕТ') or (dfline['T478K'] == 'НЕТ') or (dfline['P681R'] == 'НЕТ'):
    	return 'НЕТ'
    else:
        return 'не определено'


def ismos(dfline):
    if (dfline['del ER156-158'] == 'ДА') and (dfline['E484K'] == 'ДА') and (dfline['S494P'] == 'ДА'):
    	return 'ДА'
    elif (dfline['del ER156-158'] == 'НЕТ') or (dfline['E484K'] == 'НЕТ') or (dfline['S494P'] == 'НЕТ'):
    	return 'НЕТ'
    else:
        return 'не определено'


def isnw(dfline):
    if (dfline['del CY136-144'] == 'ДА') and (dfline['E484K'] == 'ДА') and (dfline['INS23598'] == 'ДА'):
        return 'ДА'
    elif (dfline['del CY136-144'] == 'НЕТ') or (dfline['E484K'] == 'НЕТ') or (dfline['INS23598'] == 'НЕТ'):
        return 'НЕТ'
    else:
        return 'не определено'


def isimp(dfline): # checking of mutations N501Y and E484K иной важный
    return ((dfline['E484K'] == 'ДА') or (dfline['N501Y'] == 'ДА'))

	


def classfin(dfline): # проверка на принадлежность к любому из штаммов для итоговой таблицы
    ibr = isbrit(dfline) # определяем принадлежность к каждому из штаммов
    isar = issar(dfline)
    isbr = isbras(dfline)
    isin = isind(dfline)
    isms = ismos(dfline)
    isnowe = isnw(dfline)
    important = isimp(dfline)
    if ibr == 'ДА':               # affiliation/non-affiliation to at least one of strains
        return '«Альфа»'
    elif isar == 'ДА':
        return '«Бета»'
    elif isbr == 'ДА':
        return '«Гамма»'
    elif isin == 'ДА':
    	return '«Дельта»'
    elif isms == 'ДА':
    	return 'B.1.1.523'
    elif isnowe == 'ДА':
        return 'B.1.1.370.1'
    elif important:
        return 'иной важный'
    elif (ibr == 'не определено' or isar == 'не определено' or isbr == 'не определено' or isin == 'не определено' or isms == 'не определено' or isnowe == 'не определено' or important == 'не определено'): 
        return 'не определено'    # the strain is not defined, if failed to find out
    else:
        return 'иной'


def get_table(df):
    df2 = df.copy()
    df2['Strain'] = df2.apply(lambda x: classfin(x), axis=1) # adding the last column to the table
    return df2


if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='Turns a table of stain-important mutations to the set of report tables.')
    parser.add_argument('input', type=str, help='Path to excel-file with necessary mutations.')
    parser.add_argument('output', type=str, help='Path to output excel file.')
    args = parser.parse_args()
    
    df = pd.read_excel(args.input)
    table = get_table(df)
    table.to_excel(args.output, index=False) 
