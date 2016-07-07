from __future__ import absolute_import

import sys
sys.path.append('/var/www/motif-Z/')
sys.path.append('/var/www/motif-Z/website')

import os
# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'website.settings')
from django.conf import settings

from celery.task import task
from celery import current_task
from celery.result import AsyncResult
from django.core.files import File
#import models
import json
import ui.models as models
from time import sleep
import src.motif_Z as motifZ
import src_SL.motifZ_SL as motifZ_SL

from celery import Celery
app = Celery('ui', backend='amqp://', broker='amqp://guest@localhost//', include=['ui.tasks'])

# Using a string here means the worker will not have to
# pickle the object when using Windows.
app.config_from_object('django.conf:settings')
app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)



@app.task(name='various.length.task')
def run_script(file_type, peps, data_summary, background_peps, result_file_path, statistics_file_path, debug_file_path, max_length, min_length, zscore_flag, binomial_prob_threshold, threshold_fg):
    
    results_list = motifZ.runScript(file_type, peps, data_summary, background_peps, result_file_path, statistics_file_path, debug_file_path, max_length, min_length, zscore_flag, binomial_prob_threshold, threshold_fg)
    
    print "results=",  results_list
    debug_file_path = results_list[3]
    normal_results_list=results_list[0:3]
   # json_results_list = json.dumps(results_list[0:3]) 
    #print "json_results=",  json_results_list 
    debug_file_path=results_list[3]
    #results_model = models.ResultsModel(task_id=run_script.request.id, motif_table=json_results_list, pep_type='Various', debug_file=debug_file_path)
    results_model=models.ResultsModel(task_id=run_script.request.id, motif_table=normal_results_list, pep_type='Various', debug_file=debug_file_path)
    results_model.save()
           
    return results_list

# Task for peptides with the same length 
@app.task(name='same.length.task')
def run_script_SL(file_type, peps, background_peps, result_file_path, debug_file_path, pep_length,
                  zscore_flag, binomial_prob_threshold, threshold_fg):
    
    result_table = motifZ_SL.runScript(file_type, peps, background_peps, result_file_path, debug_file_path, pep_length,
                                       zscore_flag, binomial_prob_threshold, threshold_fg)
    
    print "result_table",  result_table   
    print "task_id=", run_script_SL.request.id
    results_model = models.ResultsModel(task_id=run_script_SL.request.id, motif_table=result_table, pep_type='Same')
    results_model.save()
           
    return result_table
