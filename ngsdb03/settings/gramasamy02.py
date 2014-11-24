# settings/gramasamy02.py
import os
from .base import *

SECRET_KEY = os.environ["DJANGO_SECRET_KEY"]

DEBUG = True
TEMPLATE_DEBUG = DEBUG
EMAIL_HOST = "localhost"
EMAIL_PORT = 1025


DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2', # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'ngsdb03aa',                      # Or path to database file if using sqlite3.
        # The following settings are not used with sqlite3:
        'USER': 'ngsdb03',
        'PASSWORD': 'ngsdb03',
        'HOST': 'ngsdb',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
        'PORT': '',                      # Set to empty string for default.
    }
}

INSTALLED_APPS += ("debug_toolbar", "gunicorn")
INTERNAL_IPS = ("127.0.0.1",)
MIDDLEWARE_CLASSES += ("debug_toolbar.middleware.DebugToolbarMiddleware", )

DEBUG_TOOLBAR_CONFIG = {
    "INTERCEPT_REDIRECTS" : False
}


