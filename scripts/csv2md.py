import os
import logging
import argparse
import numpy as np
from collections import Counter, namedtuple

import pandas as pd
import inflect

try:
    import sys
    reload(sys)
    sys.setdefaultencoding('utf-8')
except:
    print('Was not able to change sys encoding to utf-8, probably b/c you\'re on Python 3.')
    pass

# One thing that may still need some manual work is the video thumbnail image insertion. I tried using this website to convert YouTube link to markdown http://embedyoutube.org/, however the thumbnail image has some randomness and may not always locate the title image. For the time being I took the screenshots manually...But otherwise the script should automate bulk of the content which is great! 

Publication = namedtuple('Publication', 'title authors keywords contactmail paperlink bloglink videolink')

_FUZZY_CATEGRORIES = [
    "Title",
    "Authors (full name, comma separated)",
    "key words",
    "Point of contact email address",
    "Link to paper",
    "Link to blog post (if any)",
    "Link to public video (e.g. YouTube)",
]

def get_info(data_row):
    return Publication(*[data_row[cat] for cat in _FUZZY_CATEGRORIES])

def format_pub_in_md(pub: Publication):
    pub_in_md = f'#### {pub.title}' 
    pub_in_md += f'\n[paper]({pub.paperlink})'
    if isinstance(pub.bloglink, str):
        pub_in_md += f' \| [blog post]({pub.bloglink})'
    pub_in_md += f'\n<br>  {pub.authors} \| _contact: {pub.contactmail}_'
    pub_in_md += f'<br>_keywords: {pub.keywords.lower()}_'
    # pub_in_md += f'<br>[![Video]({pub.videolink})]({pub.videolink})'
    return pub_in_md


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('--template_file', '-tf', type=str, default='digest_template.md')
    parser.add_argument('--input_csv', '-i', type=str, required=True)
    parser.add_argument('--output_md', '-o', type=str, required=True)
    parser.add_argument('--force_overwrite', '-f', action='store_true')
    args = parser.parse_args()

    n = args.digest_number
    p = inflect.engine()
    n_english = p.number_to_words(p.ordinal(n))
    logging.info('Parsing for the {} digest'.format(n_english))

    logging.info('Will save result to {}'.format(args.output_md))
    if os.path.isfile(args.output_md):
        if not args.force_overwrite:
            raise ValueError('Cannot overwrite existing output file!')

    logging.info('Loading template from {}'.format(args.template_file))
    with open(args.template_file, 'r') as f:
        md_template = f.read()

    logging.info('Reading {}'.format(args.input_csv))
    csv = pd.read_csv(args.input_csv)

    logging.info('Populating content...')
    
    # Construct a blurb for each pub in markdown
    md_blurbs_per_pub = []
    pubs = []

    for row_num, row in csv.iterrows():
        publication = get_info(row)
        pubs.append(publication)
        md_blurbs_per_pub.append(format_pub_in_md(publication))

    md_blurbs_per_pub = sorted(md_blurbs_per_pub)

    content = "\n".join(md_blurbs_per_pub)
    
    md = md_template.replace('$content$', content)

    logging.info('Saving digest markdown...')
    with open(args.output_md, 'w') as f:
        f.write(md)

    logging.info('Done!')
