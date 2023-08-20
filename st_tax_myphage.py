# Created: Mon Mar  7 20:05:46 2022
# Last changed: Time-stamp: <Last changed 2023-08-20 21:33:25 by Thomas Sicheritz, thomas>

import os, sys
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import streamlit as st

##########################
# import our own hacks
############
sys.path.insert(0, './modules')
sys.path.insert(0, './utils')
app_home = Path(__file__).resolve().parent
sys.path.insert(0, (app_home / 'utils').as_posix())
from our_st_tools import green_text, remove_made_with_streamlit, citation, save_file, our_exception_handler
    

############################################################
# web server details
UPLOAD_DIR = './uploads'
TITLE = "tax_MyPhage"
SUBTITLE = "- predicting the taxonomy of a phage"
DESCRIPTION = '⬅ Upload a phage genome in fasta format in the box to the left'
CITATION = """tax_MyPhage [...]
"""
CONTACT = """Contact: <a href="mailto:adm39@leicester.ac.uk">Andrew Millard </a> and <a href="mailto:thomassp@sund.ku.dk">Thomas Sicheritz-Pontén</a>"""


##########################
# Main page
############
st.set_page_config(page_title=TITLE, layout='wide')
remove_made_with_streamlit()
#st.markdown(f"{TITLE}:{SUBTITLE}")
st.title(TITLE)
st.subheader(SUBTITLE)
st.markdown(CONTACT, unsafe_allow_html=True)
green_text(DESCRIPTION)
example_button = st.button("Or have a look at an example run")
result_container = st.container()

##########################
# Sidebar
############
with st.sidebar:
    st.title('PhageCompass')
    form = st.form(key='my-form')
    col1, col2 = st.columns(2)
    mash_dist = col1.slider(min_value=0.1, max_value=0.5, label="distance threshold", value=0.3, step=0.05)
    
    #ann = form.text_input('Filter for annotation')
    #threshold = st.slider(min_value=0.01, max_value=0.25, label="distance threshold", value=0.15)
    uploaded_file = form.file_uploader('')
    #uploaded_file = form.file_uploader('Upload a phage genome in (gzipped) fasta format')
    submit = form.form_submit_button('Submit')
    citation(CITATION)


##########################
# What to run 
############
#@our_exception_handler
def run(file, result_container, verbose=False):
    result_container.empty()
    from tax_myPHAGE_nonviridic_web import Run_tax_myphage
    threads = 1
    HOME = os.path.dirname(__file__)
    base = os.path.basename(file).removesuffix('.gz').removesuffix('.fasta').removesuffix('.fna').removesuffix('.fsa')
    prefix = ''
    verbose = 1

    
    mash_df, im = Run_tax_myphage(file, base, prefix, mash_dist, threads, HOME, verbose=verbose)
    st.dataframe(mash_df)
    st.pyplot(im)
    

###############################
# wait for user clicking submit
############
if submit:
    localfile = save_file(uploaded_file, UPLOAD_DIR)
    run(localfile, result_container)

if example_button:
    file = app_home / 'GDA29J.fa'
    run(file.as_posix(), result_container)
    
