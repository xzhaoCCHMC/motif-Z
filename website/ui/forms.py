from django import forms
from django.shortcuts import redirect
from django.views.generic import FormView
from django.forms.widgets import TextInput
from django.forms.widgets import RadioSelect
from django.forms.widgets import PasswordInput
from django.utils.safestring import mark_safe

#Use customized FileField to validate size and file type
from fileUploadCheck import RestrictedFileField

class HorizRadioRenderer(RadioSelect.renderer):
    """ this overrides widget method to put radio buttons horizontally
        instead of vertically.
    """
    def render(self):
            """Outputs radios"""
            return mark_safe(u'\n'.join([u'%s\n' % w for w in self]))

class HorizRadioSelect(RadioSelect):
    renderer = HorizRadioRenderer

# Login form for the login page
class LoginForm(forms.Form):
    
    username = forms.CharField(widget=TextInput, max_length=255)
    password = forms.CharField(widget=PasswordInput)
    

class UploadFileForm(forms.Form):
    
    FILE_FORMAT_CHOICES = (
        ('csv', 'csv'),
        ('txt', 'txt'),
                                            )
    BACKGROUND_CHOICES = (
	('Uniprot', 'Uniprot human peptides'),
	('Random', 'Random peptides'),
	('UseOwn', 'Upload background from local file'),
			)

   # inputFile = forms.FileField(label='Foreground', help_text='max. 10 megabytes', required=False)
    inputFile = RestrictedFileField(content_types = ['text/plain', 'text/csv'], max_upload_size = 5242880, required = False)
    backgroundChoice = forms.ChoiceField(label='Choose background', initial='Uniprot', choices=BACKGROUND_CHOICES, required=True)
   # backgroundFile = forms.FileField(label='Background', help_text='max. 10 megabyte', required=False)
    backgroundFile = RestrictedFileField(content_types = ['text/plain'], max_upload_size = 11242880, required = False)
    fileType = forms.ChoiceField(label='File format', widget=HorizRadioSelect, choices=FILE_FORMAT_CHOICES, required=False)
    inputText = forms.CharField(label='Paste data in text box (max. 4000 peptides)', help_text='using text format',
                                widget=forms.Textarea(attrs={'row':10000, 'col':25}), required=False)
    
class UploadForegroundFileForm_SL(forms.Form):

    FILE_FORMAT_CHOICES = (
        ('csv', 'csv'),
        ('txt', 'txt'),
                                            )
    inputFile = RestrictedFileField(content_types = ['text/plain', 'text/csv'], max_upload_size = 5242880, required = False)
    fileType = forms.ChoiceField(label='File format', widget=HorizRadioSelect, choices=FILE_FORMAT_CHOICES, required=False)
    inputText = forms.CharField(label='Paste data in text box (max. 4000 peptides)', help_text='using text format',
                                widget=forms.Textarea(attrs={'row':10000, 'col':25}), required=False)

class UploadBackgroundFileForm_SL(forms.Form):
    
    SL_BACKGROUND_CHOICES = (
	('Uniprot', 'Uniprot human peptides'),
        ('Random', 'Random peptides'),
        ('UseOwn', 'Upload background from local file'),
                        )

    backgroundChoice = forms.ChoiceField(label='Choose background', initial='Uniprot', choices=SL_BACKGROUND_CHOICES, required=True)
    backgroundFile_SL = RestrictedFileField(content_types = ['text/plain'], max_upload_size = 10242880, required = False)


    
class ParametersForm(forms.Form):

    LENGTH_CHOICES = (
        ('7', '7'),
        ('8', '8'),
        ('9', '9'),
        ('10', '10'),
        ('11', '11'),
        ('12', '12'),
        ('13', '13'),
        ('14', '14'),
                     )
    Z_SCORE_CHOICES = (
                       ('1', '1'),
                       ('2', '2'),
                       ('3', '3'),
                       ('4', '4'),
                       ('5', '5'),
                       ('6', '6'),
                       ('0', 'Not use Z-Score')
                                    )

    maxLength = forms.ChoiceField(label='max. peptide length', initial='14', choices=LENGTH_CHOICES, required=True)
    minLength = forms.ChoiceField(label='min. peptide length', initial='7', choices=LENGTH_CHOICES, required=True)
    zScore = forms.ChoiceField(label='Threshold of Z-Score', initial=4, choices=Z_SCORE_CHOICES, required=True)
    bp = forms.FloatField(label='Threshold of binomial probability (default using Z-Score)', required=False)
    nf = forms.IntegerField(label='Minimum occurrence', initial=20, required=True)
    jobname = forms.CharField(label='Jobname', max_length=25, required=False)
    email = forms.EmailField(label='If you want to notified by email, input your email address (optional):', required=False)
    
class ParametersForm_SL(forms.Form):

    LENGTH_CHOICES = (
        ('7', '7'),
        ('8', '8'),
        ('9', '9'),
        ('10', '10'),
        ('11', '11'),
        ('12', '12'),
                     )
    Z_SCORE_CHOICES = (
                       ('1', '1'),
                       ('2', '2'),
                       ('3', '3'),
                       ('4', '4'),
                       ('5', '5'),
                       ('6', '6'),
                       ('0', 'Not use Z-Score')
                                    )

    pepLength = forms.ChoiceField(label='Peptide length', initial='9', choices=LENGTH_CHOICES, required=True)
    zScore = forms.ChoiceField(label='Threshold of Z-Score', initial=4, choices=Z_SCORE_CHOICES, required=True)
    #bp = forms.FloatField(label='Threshold of binomial probability (default using Z-Score)', widget=forms.HiddenInput(), initial=1000, required=False)
    bp = forms.FloatField(label='Threshold of binomial probability (default using Z-Score)', required=False)
    nf = forms.IntegerField(label='Minimum occurrence', initial=20, required=True)
    jobname = forms.CharField(label='Jobname', max_length=25, required=False)
    email = forms.EmailField(label='If you want to notified by email, input your email address (optional):', required=False)
