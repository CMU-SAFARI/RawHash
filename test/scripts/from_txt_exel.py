import pandas as pd
import openpyxl as open
import sys
import xlsxwriter as xls

file = sys.argv[1]
#n = sys.argv[2]
e = sys.argv[2]
#e = sys.argv[3]

df1 = pd.read_fwf(file, header=None)
#df = pd.DataFrame({'Col1': {df1}})
df1 = df1.replace(',','', regex=True)
#for e in {3:7}: 
 #   for n in {2:5}:
        #df = pd.read_fwf("python_eval_n{2}_e{4}.txt")
#writer = pd.ExcelWriter('RawBlend_Evaluation_bl_reg_'+str(n)+'_'+str(e)+'.xlsx', engine='xlsxwriter')
writer = pd.ExcelWriter('RawBlend_Evaluation_M1_'+str(e)+'.xlsx', engine='xlsxwriter')
df1.T.to_excel(writer, sheet_name='TestSheet', index=False)
workbook  = writer.book
worksheet = writer.sheets['TestSheet']
format1 = workbook.add_format({'num_format': '0.00000000000000000'})
worksheet.set_column(0,19, None, format1)  # Adds formatting to column C
writer.save()