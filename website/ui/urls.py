from django.conf.urls import patterns, url
from django.conf.urls import include

from django.contrib import admin
admin.autodiscover()

from ui import views

urlpatterns = patterns('',
    # ex: /polls/
    url(r'^home/$', views.home_page, name='home'),
    url(r'^vl$', views.home_variousLength, name='home_vlen'),
    url(r'^sl$', views.home_sameLength, name='home_slen'),
    url(r'^info_maxLength$', views.info_maxLength, name='info_maxLength'),
    url(r'^info_minLength$', views.info_minLength, name='info_minLength'),
    url(r'^info_pep_paste$', views.info_pep_paste, name='info_pep_paste'),
    url(r'^info_zScore$', views.info_zScore, name='info_zScore'),
    url(r'^info_bp$', views.info_bp, name='info_bp'),
    url(r'^info_nf$', views.info_nf, name='info_nf'),
    url(r'^info_pepLength$', views.info_pepLength, name='info_pepLength'),
    url(r'^info_foregroundFormat$', views.info_foregroundFormat, name='info_foregroundFormat'),
    url(r'^info_backgroundFormat$', views.info_backgroundFormat, name='info_backgroundFormat'),
    url(r'^faq/$', views.faq, name='faq'),
    url(r'^base/$', views.base, name='base'),
    url(r'^reference/$', views.reference, name='reference'),
    url(r'^algorithm/$', views.algorithm, name='algorithm'),
    url(r'^faq/$', views.faq, name='faq'),
    url(r'^about/$', views.about, name='about'),
    url(r'^poll_state/$', 'ui.views.poll_state', name='poll_state'),
    url(r'^results/task_id=(?P<task_id>.+)/jobname=(?P<jobname>.+)$', views.resultsLink, name='results'),
    url(r'^VLtask/task_id=(?P<task_id>.+)$', views.submitted_vl, name='submitted_vl'), 
    url(r'^SLtask/task_id=(?P<task_id>.+)$', views.submitted_sl, name='submitted_sl'),
    
    url(r'^login/$', views.login_page, name='login'),
    url(r'^logout/$', views.logout_page, name='logout'),
       
#     url(r'^static/(?P<path>.*)$', 'django.views.static.serve',
#         {'document_root': ''}),
#     url(r'^$', views.home, name='home'),
   

)

