#coding: utf-8

from django import forms
from django.shortcuts import redirect, render, render_to_response
from django.views.generic import FormView
from django.template import RequestContext, loader
from django.core.files.base import ContentFile
from django.http import HttpResponseRedirect, HttpResponse
from django.views.decorators.csrf import csrf_exempt
from django.core.mail import send_mail
from django.contrib.auth import authenticate, login
import json

from celery.result import AsyncResult
import tasks
from django.core.urlresolvers import reverse
#from django.utils import simplejson as json

import sys, os
import time
import decimal

sys.path.append('/var/www/motif-Z/')
import src.readData_from_UI as rd
import src_SL.readFile_SL as SL_rf
import forms as f
import models as m
from models import ResultsModel

# For upload the file and show file path
#from ajaxuploader.views import AjaxFileUploader


# for debug
import pdb

#Send email without change settings.py
import smtplib
from email.mime.text import MIMEText as text
def send_email(subject, message, sender,receiver):
    sender = sender
    receivers = receiver
    m = text(message)
    m['Subject'] = subject
    m['From'] = sender
    m['To'] = receiver

    try:
        smtpObj = smtplib.SMTP('localhost')
        smtpObj.sendmail(sender, receivers, str(m))
        print "Successfully sent email"
    except SMTPException:
        print "Error: unable to send email" 

@csrf_exempt
def poll_state(request):
	data = 'Fail'

	if request.is_ajax():
		if 'task_id' in request.POST.keys() and request.POST['task_id']:
			task_id = request.POST['task_id']
			email = request.POST['email']
			jobname = request.POST['jobname']
#			print "email==", email
			task = AsyncResult(task_id)
			print "task.state=", task.state
			if task.state == 'STARTED':
				task_state = 'Running'
				data = 'Running'
			elif task.state == 'PENDING' or task.state == 'RETRY':
				task_state = 'Waiting'
				data = 'Pending'
			elif task.state == 'SUCCESS':
				task_state = 'Finished'
				if task.result:
					data = task.result
				else:
					data = 'None'
				if email:
					send_email('Your submission at Motif-Z finished', 'Your job '+jobname+' at Motif-Z has finished, please check your result using the following url:'+'https://research.cchmc.org/motifZ/results/task_id='+task_id+'/jobname='+jobname+'\n\nIf you have any questions about motif-Z results,DO NOT reply this email, contact us by email provided on the website.', 'motifZ-admin@research.cchmc.org', email)
				else:
					print "not email received"
			else:
				task_state = task.state
				data = 'Error'
				print 'data status =', task_state
		else:
			task_state = task.state
			data = 'Error'
	else:
		task_state = task.state
		data = "Error"
		
	json_data = json.dumps({'task_state':task_state, 'task_data':data})
	print 'data check=', task_state
	return HttpResponse(json_data, mimetype='application/json')

# Login page
def login_page(request):
	
	if request.method == 'POST':
			
		loginform=f.LoginForm(request.POST)
		request_link=request.session['request_link']
		
		if loginform.is_valid():
			username=loginform.cleaned_data['username']
			password=loginform.cleaned_data['password']
			user=authenticate(username=username, password=password)
	
			if user is not None:
				if user.is_active:
					login(request, user)
					if request_link == 'various_length':
						return redirect('ui.views.home_variousLength',)
					elif request_link == 'same_length':
						return redirect('ui.views.home_sameLength',)
					else:
						return render(request, 'ui/login_page.html', {'loginform': loginform})
			
		else:
			loginform=f.LoginForm()	
	else:
		loginform=f.LoginForm()

	return render(request, 'ui/login_page.html', {'loginform': loginform})

# Logout view
from django.contrib.auth import logout

def logout_page(request):
	logout(request)
	return render(request, 'ui/logout_page.html')

def home_page(request):
	# Without login code
	if request.method == 'POST':
		if 'various_length' in request.POST:
			return redirect('ui.views.home_variousLength',)
		if 'same_length' in request.POST:
			return redirect('ui.views.home_sameLength',)
#	if request.method == 'POST':
#		if 'various_length' in request.POST:
#			request.session['request_link'] = 'various_length'
#		elif 'same_length' in request.POST:
#			request.session['request_link'] = 'same_length'
#		else:
#			request.session['request_link'] = 'Error'
#		
#		print 'request_link =', request.session['request_link']
#		
#		return redirect('ui.views.login_page')
	
	return render(request, 'ui/home.html', )


def home_variousLength(request):

	# login block
#	if not request.user.is_authenticated():
#		return render(request, 'ui/login_error.html')
	
	success_flag = 0
	print "success_flag_0=", success_flag
	
	example_file = open('/var/www/motif-Z/input/test.txt', 'r')
	example_content = example_file.read()
	example_raw = repr(example_content)
	example_text = example_raw[1:-1]
	
	if request.method == 'POST':
		
		fileform=f.UploadFileForm(request.POST, request.FILES)
		
		if fileform.is_valid():
			if fileform.cleaned_data['fileType']:
				file_type = fileform.cleaned_data['fileType']
			else:
				file_type = 'txt'
			
			if request.FILES.get('inputFile') and request.POST.get('inputText'):
				fileform = f.UploadFileForm()
				parametersform=f.ParametersForm()
				return render(request, 'ui/page_VL_file_error.html')
			
			elif request.FILES.get('inputFile'):
		# File from request will be saved on disk in this step. 
				pep_file = m.InputFile(inputfile = request.FILES['inputFile'])
		#The file object will be saved to the location specified by the upload_to argument of the corresponding FileField in model. 
				pep_file.save() 
		#Any InputFile instance has a inputfile attribute that you can use to get at the details of the uploaded file
				input_pep_file_path = pep_file.inputfile.path
				 
				data_flag = 'FILE'
				
			elif request.POST.get('inputText'):
				pep_input = request.POST['inputText']
				data_flag = 'PASTE'
		
			else:
				fileform = f.UploadFileForm()
				parametersform=f.ParametersForm()
				return render(request, 'ui/page_VL_file_error.html')
 
			
			background_type = fileform.cleaned_data['backgroundChoice']
			background_file_path = '/var/www/motif-Z/bg_files/BG_1M_VL.txt'
			if background_type == 'UseOwn' and request.FILES.get('backgroundFile'):
				background_file = m.BackgroundInputFile(backgroundFile = request.FILES['backgroundFile'])
				background_file.save()
				background_file_path = background_file.backgroundFile.path
				print "background_file_path =", background_file_path
			elif background_type != 'UseOwn' and request.FILES.get('backgroundFile'):
				#render(request, 'ui/page_VL_file_error.html')
				return HttpResponse('<h1>If upload backround from local file, choose the correct option and try again.</h1>')
			elif background_type == 'Random':
				background_file_path = '/var/www/motif-Z/bg_files/1M_randomMHC_7to14mers_wc4.txt'
			elif background_type == 'Uniprot':
				background_file_path = '/var/www/motif-Z/bg_files/BG_1M_VL.txt'
			else:
				return HttpResponse('<h1>Background file uploading error, please try again.</h1>')

			success_flag += 1
			print "success_flag_1=", success_flag
		else:
			fileform = f.UploadFileForm()
		
		parametersform = f.ParametersForm(request.POST)
		
		if parametersform.is_valid():	
			max_length = parametersform.cleaned_data['maxLength']
			min_length = parametersform.cleaned_data['minLength']
			zscore_flag = parametersform.cleaned_data['zScore']
			binomial_prob_threshold = parametersform.cleaned_data['bp']
			threshold_fg = parametersform.cleaned_data['nf']
			jobname = parametersform.cleaned_data['jobname']
			email = parametersform.cleaned_data['email']

			if int(min_length) >= int(max_length):
				fileform = f.UploadFileForm()
				parametersform=f.ParametersForm()
				return render(request, 'ui/page_VL_file_error.html')
			else:
				success_flag += 1
				print "success_flag_2=", success_flag
		else:
			parametersform = f.ParrametersForm()   
		 
		if success_flag == 2:
#			 cmd = ["python", "/home/xueheng/Documents/motif-Z/src_unix/motifZ_VL.py", "--inputf", media_pep_file_path]
#			 process_motifZ = subprocess.Popen(cmd,
#											   stderr=subprocess.STDOUT,
#											   stdout=subprocess.PIPE,
#											   shell=True)
			
			RESULT_FILE_FOLDER = '/var/www/motif-Z/results' + time.strftime('/%Y/%m/%d/')
			if not os.path.exists(RESULT_FILE_FOLDER):
				os.makedirs(RESULT_FILE_FOLDER)
				os.chmod(RESULT_FILE_FOLDER, 0777) 
			result_file_path = RESULT_FILE_FOLDER + 'file_motif_' + time.strftime("%H%M%S.txt")
			statistics_file_path = RESULT_FILE_FOLDER + 'file_stat' + time.strftime("%H%M%S.txt")
			debug_file_path = RESULT_FILE_FOLDER + 'file_debug' + time.strftime("%H%M%S.txt")
		
			# Read data from input file
			if data_flag == 'FILE':
				peps, data_summary = rd.readFile_from_UI(file_type, input_pep_file_path)
			else:
				peps, data_summary = rd.readInput_from_UI(pep_input)
			
			# Read into background file
			background_peps, bg_summary = rd.readFile_from_UI('txt', background_file_path)
			
			print "length_bg=", len(background_peps)
			task_id = 'Fail'
			job = tasks.run_script.delay(file_type, peps, data_summary, background_peps, result_file_path, statistics_file_path, debug_file_path, max_length, min_length, zscore_flag, binomial_prob_threshold, threshold_fg)
			
			status_url = 'VLtask/task_id=' + job.id
			request.session['email'] = email
			request.session['jobname'] = jobname
			return redirect(status_url)
		
	else:
		fileform = f.UploadFileForm()
		parametersform=f.ParametersForm()
		
	return render(request, 'ui/page_VL.html', {'fileform': fileform, 'parametersform': parametersform, 'exampletext': example_text})

	
def home_sameLength(request):

	# login block
#	if not request.user.is_authenticated():
#		return render(request, 'ui/login_error.html')

	success_flag = 0
	print "success_flag_0=", success_flag
	
	example_file = open('/var/www/motif-Z/input/example_SL.txt', 'r')
	example_content = example_file.read()
	example_raw = repr(example_content)
	example_text = example_raw[1:-1]
	
	if request.method == 'POST':
		
		fileform_SL=f.UploadForegroundFileForm_SL(request.POST, request.FILES)
		
		if fileform_SL.is_valid():
			file_type = fileform_SL.cleaned_data['fileType']
			
			if request.FILES.get('inputFile') and request.POST.get('inputText'):
				fileform_SL = f.UploadForegroundFileForm_SL()
				background_fileform_SL = f.UploadBackgroundFileForm_SL()
				parametersform_SL=f.ParametersForm_SL()
				return render(request, 'ui/page_SL_file_error.html')
			
			elif request.FILES.get('inputFile'):
				pep_file = m.InputFile(inputfile = request.FILES['inputFile']) 
				pep_file.save() 
				input_pep_file_path = pep_file.inputfile.path
				#print "pep_file_path =", media_pep_file_path
				data_flag = 'FILE'
				
			elif request.POST.get('inputText'):
				pep_input = request.POST['inputText']
				data_flag = 'PASTE'
		
			else:
				fileform_SL = f.UploadForegroundFileForm_SL()
				background_fileform_SL = f.UploadBackgroundFileForm_SL()
				parametersform_SL=f.ParametersForm_SL()
				return render(request, 'ui/page_SL_file_error.html')
 
			success_flag += 1
			print "success_flag_1=", success_flag
		else:
			fileform_SL = f.UploadForegroundFileForm_SL()
			
		parametersform_SL = f.ParametersForm_SL(request.POST)
		
		if parametersform_SL.is_valid():	
			pep_length = parametersform_SL.cleaned_data['pepLength']
			zscore_flag = parametersform_SL.cleaned_data['zScore']
			binomial_prob_threshold = parametersform_SL.cleaned_data['bp']
			threshold_fg = parametersform_SL.cleaned_data['nf']
			email = parametersform_SL.cleaned_data['email']
			jobname = parametersform_SL.cleaned_data['jobname']
			success_flag += 1
			print "success_flag_2=", success_flag
		else:
			parametersform_SL = f.ParametersForm_SL()   
	
		backgroundfileform=f.UploadBackgroundFileForm_SL(request.POST, request.FILES)
                if backgroundfileform.is_valid():
			background_type = backgroundfileform.cleaned_data['backgroundChoice']
			print "bk_tpye=", background_type
                        background_file_path = ''
                        if background_type == 'UseOwn' and request.FILES.get('backgroundFile_SL'):
                                background_file = m.BackgroundInputFile(backgroundFile = request.FILES['backgroundFile_SL'])
                                background_file.save()
                                background_file_path = background_file.backgroundFile.path
                                print "background_file_path =", background_file_path
                        elif background_type != 'UseOwn' and request.FILES.get('backgroundFile_SL'):
                                #render(request, 'ui/page_VL_file_error.html')
                                return HttpResponse('<h1>If upload backround from local file, choose the correct option and try again.</h1>')
                        elif background_type == 'Random' or background_type == 'Uniprot':
				if pep_length == '9':
					background_file_path = '/var/www/motif-Z/bg_files/pep_bg_9mers_1Million.txt'
				elif pep_length == '10':
                                        background_file_path = '/var/www/motif-Z/bg_files/pep_bg_10mers_100K.txt'
				elif pep_length == '11':
                                        background_file_path = '/var/www/motif-Z/bg_files/pep_bg_11mers_100K.txt'
				elif pep_length == '12':
                                        background_file_path = '/var/www/motif-Z/bg_files/pep_bg_12mers_100K.txt'
				elif pep_length == '7':
                                        background_file_path = '/var/www/motif-Z/bg_files/pep_bg_flank_7mers_100K.txt'
				elif pep_length == '8':
                                        background_file_path = '/var/www/motif-Z/bg_files/pep_bg_8mers_100K.txt'
				else:
					return HttpResponse('<h1>Unrecoverable error in uploading background file.</h1>')
                        else:
                                return HttpResponse('<h1>Background file uploading error, please try again.</h1>')
		else:
			backgroundfileform = f.UploadBackgroundFileForm_SL()
	
		if success_flag == 2:

			RESULT_FILE_FOLDER = '/var/www/motif-Z/results' + time.strftime('/%Y/%m/%d/')
			if not os.path.exists(RESULT_FILE_FOLDER):
				os.makedirs(RESULT_FILE_FOLDER) 
				os.chmod(RESULT_FILE_FOLDER, 0777)
			result_file_path = RESULT_FILE_FOLDER + 'file_motif_' + time.strftime("%H%M%S.txt")
			debug_file_path = RESULT_FILE_FOLDER + 'file_debug' + time.strftime("%H%M%S.txt")
			
			# Read data from input file
			if data_flag == 'FILE':
				peps, data_summary = rd.readFile_from_UI(file_type, input_pep_file_path)
				print "input_pep_file_path=", input_pep_file_path
			else:
				peps, data_summary = rd.readInput_from_UI(pep_input)
		 	# Read into background file
                        background_peps = SL_rf.read_txtPep(background_file_path)
	
			task_id = 'Fail'
			job = tasks.run_script_SL.delay(file_type, peps, background_peps,  result_file_path, debug_file_path, pep_length,
										 zscore_flag, binomial_prob_threshold, threshold_fg)
			request.session['email'] = email
			request.session['jobname'] = jobname
			status_url = 'SLtask/task_id=' + job.id
			return redirect(status_url)

	else:
		fileform_SL = f.UploadForegroundFileForm_SL()
		backgroundfileform = f.UploadBackgroundFileForm_SL()
		parametersform_SL=f.ParametersForm_SL()
		
	return render(request, 'ui/page_SL.html', {'fileform_SL': fileform_SL, 'parametersform_SL': parametersform_SL, 'backgroundfileform': backgroundfileform, 'exampletext': example_text})

# Submitted pages
def submitted_vl(request, task_id):
	email = request.session.get('email')
	jobname = request.session.get('jobname')
	return render(request, 'ui/submitted.html', {'task_id': task_id, 'email': email, 'jobname':jobname}) 
 
def submitted_sl(request, task_id):
	email = request.session.get('email')
	jobname = request.session.get('jobname')
	return render(request, 'ui/submitted_SL.html', {'task_id': task_id, 'email': email, 'jobname':jobname}) 

# Information page for parameters
def info_maxLength(request):
	return render(request, 'ui/info_maxLength.html', )
def info_minLength(request):
	return render(request, 'ui/info_minLength.html', )
def info_pep_paste(request):
	return render(request, 'ui/info_pep_paste.html', )
def info_zScore(request):
	return render(request, 'ui/info_zScore.html', )
def info_bp(request):
	return render(request, 'ui/info_bp.html', )
def info_nf(request):
	return render(request, 'ui/info_nf.html', )
def info_pepLength(request):
	return render(request, 'ui/info_pepLength.html', )
def info_foregroundFormat(request):
	return render(request, 'ui/info_foregroundFormat.html')
def info_backgroundFormat(request):
        return render(request, 'ui/info_backgroundFormat.html')


# Information page for home page
def base(request):
	return render(request, 'ui/base.html', )

def faq(request):
	return render(request, 'ui/faq.html', )

# Information page for home page
def reference(request):
	return render(request, 'ui/reference.html', )

def algorithm(request):
	return render(request, 'ui/algorithm.html', )

def about(request):
	return render(request, 'ui/about.html', )

def resultsLink(request, task_id, jobname):
	jobname = jobname
	result_data = ResultsModel.objects.get(task_id=task_id)
	task_id = result_data.task_id
	pep_type = result_data.pep_type
	motif_table = json.dumps(result_data.motif_table, default=decimal_default)
	return render(request, 'ui/results.html', {'task_id':task_id, 'pep_type':pep_type, 'motif_table':motif_table, 'jobname':jobname})

def decimal_default(obj):
    if isinstance(obj, decimal.Decimal):
        return float(obj)
    raise TypeError
