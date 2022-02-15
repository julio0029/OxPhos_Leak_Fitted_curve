#!/usr/bin/env python3
#-*- coding: utf-8 -*-

'''-------------------------------------------------------------------------------
Copyright© 2021 Jules Devaux / Alice Harford. All Rights Reserved
Open Source script under Apache License 2.0
----------------------------------------------------------------------------------

Fit respiration curve for the Leak/OxPhos study and extract relevant parameters.
'''


import os
current_path = os.path.dirname(os.path.abspath(__file__))


#======== PARAMETERS =========

CHAMBER_VOLUME = 2 #ml
WINDOW = 20 #Moving averages
MAX_kPA = 24 #kPa
GRAPHING = False
_saving = True # SAve summary file and stats
FOLDER = f'{current_path}/CSV' # Folder containing all csvs
FILENAME = None#'Leak_OXPHOS_18_Trial10.csv' #Fro debugging

#=============================


import sys
import numpy as np
import pandas as pd
import pingouin as pg
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

if GRAPHING is True:
	try: from graph import Graph
	except:
		print(f"graph.py required in current directory {current_path}")



def hill(S, Vmax, Km, n=1):
	return (Vmax*S**n)/(Km+S**n)


def fit_hill_curve(x, y):
	'''
	x_y is df with x as index, y as values
	'''
	Vmax = y.max()
	xmax, xmin = x[-1], x[0]
	dx = x[-1]-x[0]

	# parameters as [max, min, dt, ]
	param_initial = [Vmax, 0.5*dx, 1]
	param_bounds = ([Vmax-0.5*Vmax, xmin*0.1, 0.01],
					[Vmax+0.5*Vmax, xmax*10, 100])
	bounds = [(Vmax-0.5*Vmax, Vmax+0.5*Vmax),
			(xmin, xmax),
			(0.01, 100)]


	popt, _ = curve_fit(
		hill, x, y,
		p0 = param_initial,
		bounds=param_bounds)

	Vmax, pP50, _hill =popt[0], popt[1], popt[2]
	y_fit = hill(x, Vmax, pP50, _hill)

	# Goodness of fit:
	# residual sum of squares
	ss_res = np.sum((y - y_fit) ** 2)
	ss_tot = np.sum((y - np.mean(y)) ** 2)
	r2 = 1 - (ss_res / ss_tot)

	predicted=pd.Series(y_fit, name='predicted')

	return Vmax, pP50, _hill, predicted, r2


def get_JO2(PO2, mass=5, _resample=0.1):
	'''
	Takes a csv file and pre-process to fit curve
	Requires mass !
	- smoothes PO2
	- recalculate JO2
	'''

	#index needs to be timed by 2 since recording time=2s
	_time=PO2.index*2

	PO2=pd.DataFrame(PO2.values, index=_time, columns=['PO2'])

	# Moving average
	PO2=PO2.rolling(window=int(WINDOW/2)).mean()

	# Calculate JO2
	JO2=(((PO2.shift()-PO2)*CHAMBER_VOLUME)/mass)*100
	JO2=pd.DataFrame(JO2.values, index=_time,columns=['JO2'])

	# Create df with PO2 as index
	JO2['PO2']=PO2
	JO2=JO2.loc[JO2['PO2']>=0]
	JO2=JO2.loc[JO2['JO2']>=0]
	JO2=JO2.drop_duplicates(subset='PO2',keep='first')

	JO2 = JO2.set_index('PO2'
			).sort_index(
			).rename(columns={JO2.columns[0]:'JO2'}
			).replace([np.inf, -np.inf, np.nan, float('nan')], 0, inplace=False
			).rolling(window=WINDOW).mean(
			).dropna(
			).astype('float64')

	return JO2


def extract_csvs():
	# All chambers saved as dict
	chambers = []
	for temperature in [18, 20, 24, 25, 26, 27, 30]:
		_path = f"{FOLDER}/{temperature}/"
		for file in os.listdir(_path):

			if '.csv' in file:
				print(f'Exctracting {file}...')
				df=pd.read_csv(f"{_path}{file}", low_memory=False)
				
				
				df['Event Name']=df['Event Name'].fillna(method='ffill')

				for chb in ['A:', 'B:']:
					col=[c for c in df.columns if chb in c][0]
					idx_col = df.columns.get_loc(col)
					mass = df.iloc[0,idx_col+1]

					# Select only to Anoxia part
					# Some files had the whole experiments
					if 'S' in df['Event Name'].values:
						select_df=df[df['Event Name']=='S']
						#print(f"{file} with S")
					else:
						select_df=df
						#print(f"{file} without S")
					
					select_df=select_df.loc[select_df[col]<=MAX_kPA].sort_index()
					PO2=select_df.loc[:,col]

					# Exctract JO2
					kPa_to_nM=(224.47/20.37)
					JO2=get_JO2(PO2, mass)*kPa_to_nM*30

					# Define OXPHOS - LEAK
					if chb[0] == 'A':
						_state = 'OXPHOS'
					else: _state = 'LEAK'


					# Store chamber parameters into chambers list
					chambers.append({
						'filename':file,
						'temperature':temperature,
						'chamber':chb[0],
						'state':_state,
						'PO2':PO2,
						'PO2max':PO2.max(),
						'JO2':JO2,
						'expJO2max':JO2.max(),
						'mass':mass
						})

	return chambers


def process_all_chambers(chambers):
	_dicts=[] #temprary list to avoid overwriting chambers list
	errors=[]
	for chamber in chambers:
		print(f"Processing {chamber['filename']}")
		#try:
		if True:

			# Fit curve to experimental PO2 
			fJO2max, fP50, _hill, fitted, r2 = fit_hill_curve(chamber['JO2'].index.values, chamber['JO2']['JO2'].values)

			# Extend curve to 0 -> 100% PO2
			# Standardise PO2 for comaprison between samples
			X=np.arange(0, MAX_kPA, 0.01)
			predicted=hill(X, fJO2max, fP50, _hill)

			# plt.scatter(X, predicted, label='pred')
			# plt.scatter(chamber['JO2'].index,chamber['JO2']['JO2'].values, label='expe')
			# plt.show()

			# JO2max at 100% saturation 
			kPa21_JO2max=predicted.max()


			# Extract P50 from fitted curve ()
			JO2_50=kPa21_JO2max/2
			cP50=_hill*(((JO2_50*fP50)/(kPa21_JO2max-JO2_50))**(_hill-1))


			# Calculate area under the curve from Standardised
			auc=np.trapz(predicted, X, 0.1)


			# Update parameters for chamber
			chamber.update(
				{'JO2': chamber['JO2'],
				'fJO2max': fJO2max,
				'fP50':fP50,
				'hill':_hill,
				'kPa21_JO2max':kPa21_JO2max,
				'cP50':cP50,
				'goodness of fit':r2,
				'area under curve':auc,
				'predicted':predicted
				})

			_dicts.append(chamber)

		# except Exception as e:
		# 	print(f"{chamber['filename']}: {e}")
		# 	errors.append(chamber['filename'])
	for e in errors:print(e)
	
	return _dicts



def create_summary(chambers):

	# Check outliers
	outliers=[]
	for chamber in chambers:
		if float(chamber['goodness of fit'])<=0.5:
			outliers.append({
				'filename':chamber['filename'],
				'chamber': chamber['chamber']
				})

		print("============= OUTLIERS ==============")
		for o in outliers:
			print(f"{o['filename']} - chamber {o['chamber']}")
		print("\n\n")

		# Create summary table
		summary=[]
		for i in range(len(chambers)):
			row={k:v for k,v in chambers[i].items() if k not in ['PO2','predicted','JO2']}

			row=pd.DataFrame(row, index=[i])
			summary.append(row)
		summary=pd.concat(summary)

		# Save summary
		if _saving is True:
			summary.to_csv('summary.csv')

		return summary


def do_stats(df):
	
	# Create Excel template
	fileName=f'Summary.xlsx'
	writer=pd.ExcelWriter(fileName, engine='xlsxwriter')

	# First, add data to first tab
	df.to_excel(writer, sheet_name='main') 


	between=['temperature','state']
	parameters=['PO2max',
				'fJO2max','fP50','hill','kPa21_JO2max',
				'cP50','goodness of fit','area under curve']


	#===== Do Averages and SEM
	summary=df.groupby(between).agg(['mean','sem'])
	summary.to_excel(writer, sheet_name='averages') 


	#===== Do two-Way ANOVA ====
	aov_row=0 #Keep track of rows to append ANOVA tests
	hoc_row=0

	# Loop through parameter for ANOVA
	for parameter in parameters:
		aov=pg.anova(df, dv=parameter,
				between=between,
				detailed=True)
		aov['Parameter']=parameter

		# Append to Excel, spacing of 3 rows
		aov.to_excel(writer, sheet_name="ANOVAs", startrow=aov_row , startcol=0)   
		aov_row+=(len(aov.index)+3)

		# ==== POST-HOCs
		groups=[{'groups': 'state', 'between': 'temperature'},
				{'groups': 'temperature', 'between': 'state'}]
		for g in groups:
			for group, gdf in df.groupby(g['groups']):
				# Do the test
				ph=pg.pairwise_tukey(data=gdf, dv=parameter, between=g['between'])
				ph['Parameter']=parameter
				ph['group']=group
		
				# Append to Excel
				ph.to_excel(writer, sheet_name="PostHocs", startrow=hoc_row , startcol=0)   
				hoc_row+=(len(ph.index)+3)


	writer.save()
	print(f'Saved {fileName}')





def main():

	chambers=extract_csvs()
	chambers=process_all_chambers(chambers)
	df=create_summary(chambers)

	#df=pd.read_csv('summary.csv', index_col=0)
	do_stats(df)





	# select one file to graph
	if GRAPHING is True:
		if FILENAME is not None:
			chambers=[c for c in chambers if c['filename']==FILENAME]

		for c in chambers:	
			JO2=c['JO2']
			predicted=c['predicted']
			Graph().graph(
				[[JO2.index,
				{'y': JO2['JO2'], 'type': 'line', 'color': 'red'},
				{'y': predicted, 'type': 'line', 'color': 'blue'}]],
				title=f"{c['temperature']}-{c['chamber']}")

if __name__ == '__main__':
	main()

