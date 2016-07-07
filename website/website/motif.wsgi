import os
import sys

path = '/var/www/motif-Z/website'
if path not in sys.path:
    sys.path.append(path)
sys.path.append('/usr/lib/django2.6/site-packages')

os.environ['DJANGO_SETTINGS_MODULE'] = 'website.settings'

import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()

