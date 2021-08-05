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
WINDOW = 20
GRAPHING = False
FOLDER = f'{current_path}/CSV'

#=============================


import sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

if GRAPHING is True:
	try: from graph import Graph
	except:
		print(f"graph.py not in current directory {current_path}")
		print("Importing pyplot, yet this code needs to be updated to plot")
		from matplotlib import pyplot as plt




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

	return popt



def get_JO2(PO2, mass=5):
	'''
	Takes a csv file and pre-process to fit curve
	Requires mass !
	- smoothes PO2
	- recalculate JO2
	'''

	_time=PO2.index*2

	PO2=pd.DataFrame(PO2.values, index=_time, columns=['PO2'])
	#PO2=PO2.loc[PO2['PO2']>0]

	# Moving average
	PO2=PO2.rolling(window=WINDOW).mean()

	# Calculate JO2
	JO2=(((PO2.shift()-PO2)*CHAMBER_VOLUME)/mass)*100
	JO2=pd.DataFrame(JO2.values, index=_time,columns=['JO2'])

	# Create df with PO2 as index
	JO2['PO2']=PO2
	JO2=JO2.loc[JO2['PO2']>=0]
	JO2=JO2.loc[JO2['JO2']>=0]

	JO2 = JO2.set_index('PO2'
			).sort_index(
			).dropna(
			).rename(columns={JO2.columns[0]:'JO2'}
			).replace([np.inf, -np.inf, np.nan, float('nan')], 0, inplace=False
			).astype('float64')


	return JO2


def process_chamber(chamber):
	'''
	Chamber as {
'				'filename':str,
				'temperature'float:,
				'chamber':str,
				'PO2':df,
				'mass':float}
	'''

	JO2=get_JO2(chamber['PO2'], chamber['mass'])
	x=JO2.index.values
	y=JO2['JO2'].values

	popt = fit_hill_curve(x, y)
	Vmax, pP50, _hill =popt[0], popt[1], popt[2]

	predicted = hill(x, Vmax, pP50, _hill)
	predicted=pd.Series(predicted, name='predicted')

	return pP50, _hill, predicted, JO2



def main():

	# Create list of dictionary
	chambers = []
	for temperature in [20, 24, 25, 26, 27, 30]:
		_path = f"{FOLDER}/{temperature}/"
		for file in os.listdir(_path):

			if '.csv' in file:
				df=pd.read_csv(f"{_path}{file}")
				
				# Select the portion from S to STOP
				df['Event Name']=df['Event Name'].fillna(method='ffill')

				for chb in ['A:', 'B:']:
					col=[c for c in df.columns if chb in c][0]
					idx_col = df.columns.get_loc(col)
					mass = df.iloc[0,idx_col+1]

					select_df=df[df['Event Name']=='S']

					chambers.append({
						'filename':file,
						'temperature':temperature,
						'chamber':chb[0],
						'PO2':select_df.loc[:,col],
						'mass':mass
						})

	_dicts=[]
	for chamber in chambers:
		print(f"Processing {chamber['filename']}")
		try:
			pP50, _hill, predicted, JO2 = process_chamber(chamber)
			JO2max=predicted.max()
			# Define OXPHOS - LEAK
			if chamber['chamber'] == 'A':
				_state = 'OXPHOS'
			else: _state = 'LEAK'



			# Calculate P50 from predicted curve
			JO2_50=JO2max/2
			cP50=hill(JO2_50, JO2max, pP50, _hill)


			# Calculate area under the curve
			scope = np.trapz(predicted, JO2.index,0.1)


			# Update parameters for chamber
			chamber.update(
				{'JO2': JO2,
				'JO2max': JO2max,
				'state':_state,
				'pP50':pP50,
				'cP50':cP50,
				'hill':_hill,
				'scope':scope,
				'predicted':predicted
				})

			_dicts.append(chamber)

		except Exception as e: print(f"{chamber['filename']}: {e}")
	chambers=_dicts



	# Check outliners
	outliners=[]
	for chamber in chambers:
		if chamber['JO2max']<chamber['cP50']:
			outliners.append({
				'filename':chamber['filename'],
				'chamber': chamber['chamber']
				})



	# Create summary table
	summary=[]
	for i in range(len(_dicts)):
		row=pd.DataFrame(_dicts[i], index=[i])
		row.drop(['PO2','predicted','JO2'], axis=1, inplace=True)
		summary.append(row)
	summary=pd.concat(summary)

	# Save summary
	_saving=True
	if _saving is True:
		summary.to_csv('summary.csv')


	# Get mean and SD for each temp
	pass



	# select one file to graph
	select=None
	for c in chambers:
		if c['filename']== 'Leak_OXPHOS_20_Trial3.csv':
			select=c
			break
	JO2=select['JO2']
	predicted=select['predicted']

	GRAPHING=False

	if GRAPHING is True:
		Graph().graph(
			[[JO2.index,
			{'y': JO2['JO2'], 'type': 'line', 'color': 'red'},
			{'y': predicted, 'type': 'line', 'color': 'blue'}]])

if __name__ == '__main__':
	main()


