"""
Django settings for djangopage project.

Generated by 'django-admin startproject' using Django 3.0.7.

For more information on this file, see
https://docs.djangoproject.com/en/3.0/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/3.0/ref/settings/
"""

import os
import django_on_heroku
import dj_database_url

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import storages.backends.s3boto3
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/3.0/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = '244e123eb364ed5e4f4ffc0acc3084da43c951036847108'

# SECURITY WARNING: don't run with debug turned on in production!
#DEBUG = False
DEBUG = True
#SECRET_KEY = os.environ.get('SECRET_KEY ')
#DEBUG = (os.environ.get('DEBUG_VALUE') == "True")



ALLOWED_HOSTS = ['127.0.0.1:8000', '127.0.0.1', 'bohui119.herokuapp.com','https://database-30lc.onrender.com','127.0.0.1:5432']
# ALLOWED_HOSTS = []


# Application definition

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django_plotly_dash.apps.DjangoPlotlyDashConfig',
    'nickelcitrate.apps.NickelcitrateConfig',
    'nickelammonia.apps.NickelammoniaConfig',
    'nickelpure.apps.NickelpureConfig',
    'nickelca.apps.NickelcaConfig',
    'copper.apps.CopperConfig',
    'home.apps.HomeConfig',
    'letsencrypt',
    'storages'

    # 'crispy_forms',
    # 'dpd_static_support',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',

    #
    # 'whitenoise.middleware.WhiteNoiseMiddleware',
    #

    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'django_plotly_dash.middleware.BaseMiddleware',
    'django_plotly_dash.middleware.ExternalRedirectionMiddleware',

]

ROOT_URLCONF = 'djangopage.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': ['templates'],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'djangopage.wsgi.application'


# Database
# https://docs.djangoproject.com/en/3.0/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
    },
    'default': dj_database_url.config(
        # Feel free to alter this value to suit your needs.
        default='postgresql://postgres:postgres@127.0.0.1:8000/database', conn_max_age=600)
}

# Password validation
# https://docs.djangoproject.com/en/3.0/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/3.0/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# CRISPY_TEMPLATE_PACK = 'bootstrap4'


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/3.0/howto/static-files/

# STATICFILES_FINDERS = [
#     'django.contrib.staticfiles.finders.FileSystemFinder',
#     'django.contrib.staticfiles.finders.AppDirectoriesFinder',
#     'django_plotly_dash.finders.DashAssetFinder',
#     'django_plotly_dash.finders.DashComponentFinder',
#     'django_plotly_dash.finders.DashAppDirectoryFinder',
# ]

PLOTLY_COMPONENTS = [
    # Common components
    'dash_core_components',
    'dash_html_components',
    'dash_renderer',
    # django-plotly-dash components
    'dpd_components',
    # static support if serving local assets
    'dpd_static_support',
    # Other components, as needed
    # # 'dash_bootstrap_components',
    # 'dpd_static_support'
]

PLOTLY_DASH = {
    # Route used for the message pipe websocket connection
    "ws_route" :   "dpd/ws/channel",

    # Route used for direct http insertion of pipe messages
    "http_route" : "dpd/views",

    # Flag controlling existince of http poke endpoint
    "http_poke_enabled" : True,

    # Insert data for the demo when migrating
    "insert_demo_migrations" : False,

    # Timeout for caching of initial arguments in seconds
    "cache_timeout_initial_arguments": 60,

    # Name of view wrapping function
    "view_decorator": None,

    # Flag to control location of initial argument storage
    "cache_arguments": True,

    # Flag controlling local serving of assets
    "serve_locally": False,
}

X_FRAME_OPTIONS = 'SAMEORIGIN'

STATIC_ROOT = os.path.join(BASE_DIR, 'staticfiles')
STATIC_URL = '/static/'
# STATICFILES_LOCATION = 'static'
# STATIC_ROOT = os.path.join(BASE_DIR, 'staticfiles')

#STATICFILES_DIRS = (
#os.path.join(BASE_DIR, "static"),


django_on_heroku.settings(locals())
AWS_ACCESS_KEY_ID = os.environ.get('AWS_ACCESS_KEY_ID')
AWS_SECRET_ACCESS_KEY = os.environ.get('AWS_SECRET_ACCESS_KEY')
AWS_STORAGE_BUCKET_NAME = os.environ.get('AWS_STORAGE_BUCKET_NAME')

AWS_S3_FILE_OVERWRITE = False
AWS_DEFAULT_ACL = None

DEFAULT_FILE_STORAGE = 'storages.backends.s3boto3.S3Boto3Storage'
