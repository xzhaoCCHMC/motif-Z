from django.db import models
from json_field import JSONField
from datetime import datetime
from django.utils import timezone

# Create your models here.
from formatChecker import ContentTypeRestrictedFileField


class InputFile(models.Model):
   # inputfile = models.FileField(upload_to='inputfiles/%Y/%m/%d')
    inputfile = ContentTypeRestrictedFileField(upload_to='inputfiles/%Y/%m/%d', content_types=['text/plain', 'text/csv',],max_upload_size=1100000,blank=True, null=True)
    date_created = models.DateTimeField('date created', default=datetime.now)
    class Meta:
        app_label = 'ui'
    def __unicode__(self):
        return "InputFileModel<%s>" % (self.inputfile)
    

class BackgroundInputFile(models.Model):
   # backgroundFile = models.FileField(upload_to='bginputfiles/%Y/%m/%d')
    backgroundFile = ContentTypeRestrictedFileField(upload_to='bginputfiles/%Y/%m/%d', content_types=['text/plain',],max_upload_size=11000000,blank=True, null=True)
    date_created = models.DateTimeField('date created', default=datetime.now)
    class Meta:
        app_label = 'ui'
    def __unicode__(self):
        return "BackgroundInputFileModel<%s>" % (self.backgroundFile)

    
#Table in SQL database storing results
class ResultsModel(models.Model):

    PEPTIDE_LENGTH = (
	('V', 'Various'),
	('S', 'Same'),
    )

    task_id = models.CharField(max_length=255)
    pep_type = models.CharField(max_length=10)
    motif_table = JSONField()
    debug_file = models.CharField(max_length=255, blank=True, null=True)
    #email = models.EmailField(max_length=75, blank=True, null=True)
    date_created = models.DateTimeField('date created', default=datetime.now)
    class Meta:
	app_label = 'ui'

    def __unicode__(self):
        return "ResultsModel<%s>" % (self.task_id)
