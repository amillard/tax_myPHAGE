#!/usr/bin/env python3.9
# Created: Mon May 16 12:26:05 2022
# Last changed: Time-stamp: <Last changed 2023-08-20 20:29:53 by Thomas Sicheritz, thomas>

import os, sys
import glob, gzip
import streamlit as st

unsafe_password = '_cbdtracker_'
font='monospace'
def green_text(txt, st=st):st.markdown(f'<p style="color:#09ab3b;font-family:{font};">{txt}</p>', unsafe_allow_html=True)
def blue_text(txt, st=st):st.markdown(f'<p style="color:#0068c9;font-family:{font};">{txt}</p>', unsafe_allow_html=True)
def red_text(txt, st=st):st.markdown(f'<p style="color:#ff2b2b;font-family:{font};">{txt}</p>', unsafe_allow_html=True)
def yellow_text(txt, st=st):st.markdown(f'<p style="color:#faca2b;font-family:{font};">{txt}</p>', unsafe_allow_html=True)
def white_text(txt, st=st):st.markdown(f'<p style="color:#ffffff;font-family:{font};">{txt}</p>', unsafe_allow_html=True)
def black_text(txt, st=st):st.markdown(f'<p style="color:#000000;font-family:{font};">{txt}</p>', unsafe_allow_html=True)


def secure_filename(name):
    import unicodedata
    import re
    name = unicodedata.normalize('NFKD', str(name)).encode('ascii', 'ignore').decode('ascii')
    name = re.sub('[^a-zA-Z0-9_]','_', name)
    return name.strip('-_')


def remove_made_with_streamlit():
    # remove the made with streamlit
    hide_streamlit_style = """<style>#MainMenu {visibility: hidden;}\nfooter {visibility: hidden;}\n</style> """
    st.markdown(hide_streamlit_style, unsafe_allow_html=True)

   
def save_file(uploaded_file, UPLOAD_DIR='./uploads'):
    localfile = os.path.join(UPLOAD_DIR, uploaded_file.name)
    fid = open(localfile, 'wb+')
    fid.write(uploaded_file.getvalue())
    fid.close()
    return localfile


def citation(CITATION):
    st.markdown(f"""<details style="font-size:11px"> <summary>Citation
    <span class="icon">📖</span>
    </summary>
    <p>
    <small>{CITATION}</small>
    </p>
    </details>""", unsafe_allow_html=True)

def HELP(txt, name):
    st.markdown(f"""<details style="font-size:20px"> <summary>{name}
    <span class="icon">❓</span>
    </summary>
    <p>
    <large>{txt}</large>
    </p>
    </details>""", unsafe_allow_html=True)
    

##########################
# for background picture
############
@st.cache_data(persist=True)
def get_base64_of_bin_file(bin_file):
    import base64
    with open(bin_file, 'rb') as f:
        data = f.read()
    return base64.b64encode(data).decode()


def set_png_as_page_bg(png_file):
    bin_str = get_base64_of_bin_file(png_file)
    page_bg_img = '''
    <style>
    .stApp {
    background-image: url("data:image/png;base64,%s");
    background-size: cover;
    }
    </style>
    ''' % bin_str
    
    st.markdown(page_bg_img, unsafe_allow_html=True)
    return


def our_exception_handler(func):
    def inner_function(*args, **kwargs):
        try:
            func(*args, **kwargs)
        except:
            from PIL import Image
            image = Image.open('./assets/phage_error.png')
            
            col1, col2, col3 = st.columns(3)
            with col1:
                st.write(' ')
            with col2:
                st.image(image, width=200)
                red_text("Something went wrong ...")
            with col3:
                st.write(' ')

    return inner_function

def add_vertical_space(num_lines: int = 1):
    """Add vertical space to your Streamlit app."""
    for _ in range(num_lines):
        st.write("")

# find icons as unicode characters
# https://html-css-js.com/html/character-codes/icons/
