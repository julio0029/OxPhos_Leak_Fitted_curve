#!/usr/bin/env python3
#-*- coding: utf-8 -*-
'''-------------------------------------------------------------------------------
Copyright© 2019 Jules Devaux All Rights Reserved
This script may not be copied, altered, distributed or processed by any rights,
unless granted by the owner (i.e. Jules Devaux)
----------------------------------------------------------------------------------

'''


import datetime, time

import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


#////////////////////////////////////// GRAPH \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
class Graph():

	def __init__(self):
		self.colors={
			'reds':['#6C0020','#810026','#95002C','#AA0032','#BE0038','#D2003E','#E70044','#e91956','#f0668e','#FF0000'],
			'blues':['#309BCD','#4DA9D4','#6BB7DB','#88C5E2','#A6D4E9','#C3E2F0','#E1F0F7'],
			'greens':['#12b500','#65eaaf','#c2f1ce','#daf7e5','#a6f79d','#44f35e']}
		

	def set_line_color(self, y_label, color=None):
		if 'rel' in y_label:
			color = "#b0b0b0"
		elif "PO2" in y_label.upper(): 
			color=self.colors.get('blues')
		elif "JO2" in y_label.upper():
			color=self.colors.get('reds')
		elif 'ATP' in y_label.upper():
			color=self.colors.get('greens')

		else:
			import random
			try:color=random.choice(self.colors.get(f'{color}s'))
			except:color="#b0b0b0"

		return color


	def set_Theme(self, ax, y_label=None):
		# A few options:
		ax.grid(True)
		ax.set_facecolor('#333744')
		ax.set_autoscale_on
		ax.spines['bottom'].set_color('#808595')
		ax.spines['top'].set_color('#808595')
		ax.spines['left'].set_color('#808595')
		ax.spines['right'].set_color('#808595')
		ax.tick_params(axis='both', direction='in',
						color='#808595', labelcolor='#808595',
						grid_color='#465063', grid_linestyle='--', grid_linewidth=0.5)
		return ax


	def set_label(self,label):
		if ('JO2' in label) or ('jo2' in label):
			return 'JO2 (O2*(s*mg)-1)'
		elif ('PO2' in label) or ('po2' in label):
			return 'PO2 (kPa)'
		elif 'date' in label:
			return 'Date - Time'
		else: return 'N/A'


	def set_ax(self, ax, x_y_s):
		x=x_y_s[0]#.to_list()
		day=mdates.DayLocator()
		try:x_label=self.set_label(str(x.name))
		except:x_label='x'
		try:y_label=self.set_label(str(x_y_s[1].name))
		except:y_label=self.set_label(str(x_y_s[1]['y'].name))
		#--- sort y ---
		for j in range(1,len(x_y_s)):
			color=self.set_line_color(y_label)[j]
			if type(x_y_s[j]) is dict:
				label=x_y_s[j]['y'].name
				y=x_y_s[j]['y'].to_list()
				try: color=self.set_line_color(y_label='None', color=x_y_s[j]['color'])
				except:pass
				if x_y_s[j]['type']=='line':ax.plot(x,y,color=color,linewidth=1,label=label)
				elif x_y_s[j]['type']=='scatter':ax.scatter(x,y,color=color,label=label)
				elif x_y_s[j]['type']=='bar':ax.bar(x,y,color=color,label=label,width=0.5)
			else:ax.plot(x , x_y_s[j].to_list(), color=color, label=x_y_s[j].name)
		ax=self.set_Theme(ax)
		ax.set_xlabel(x_label, labelpad=10, color='#808595', fontsize='small',fontweight='bold')		
		ax.set_ylabel(y_label, labelpad=10, color='#808595', fontsize='x-small',fontweight='bold')
		ax.legend(loc='upper left', fontsize='xx-small')
		return ax
	

	def graph(self, x_y_s, x_range=None, title=None):
		'''
		Kwargs are passed to matplotlib plotting functions.
		'''
		if type(x_y_s[0]) is pd.core.series.Series:
			x_y_s=[x_y_s]
		subplot_nb=len(x_y_s)
		fig, axs = plt.subplots(subplot_nb, 1,facecolor='#283036', sharex=True)
		#--- one plot (need to otherwise returns error) --
		if subplot_nb==1:
			x_y_s=x_y_s[0]
			axs=self.set_ax(axs,x_y_s)
			
		#--- >1 subplot --- 
		else:
			for i in range(subplot_nb):
				axs[i]=self.set_ax(axs[i],x_y_s[i])
			plt.subplots_adjust(hspace=0)
		if x_range:
			plt.xlim(x_range)
		if title:
			plt.suptitle(title, fontsize=10, fontweight='heavy', color='#808595', alpha=0.2, x=0.5,y=0.5)
		fig.autofmt_xdate()
		plt.show()


if __name__ == '__main__':
	pass
