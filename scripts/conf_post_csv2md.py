import os
import logging
import argparse
import numpy as np
from collections import Counter, namedtuple
import csv
import inflect

try:
    import sys
    reload(sys)
    sys.setdefaultencoding('utf-8')
except:
    print('Was not able to change sys encoding to utf-8, probably b/c you\'re on Python 3.')
    pass

# One thing that may still need some manual work is the video thumbnail image insertion. I tried using this website to convert YouTube link to markdown http://embedyoutube.org/, however the thumbnail image has some randomness and may not always locate the title image. For the time being I took the screenshots manually...But otherwise the script should automate bulk of the content which is great! 

Publication = namedtuple('Publication', 'timestamp title authors abstract keywords contactmail paperlink awards bloglink videolink')

_FUZZY_CATEGRORIES = [
    "Timestamp",
    "Title",
    "Authors (full name, comma separated)",
    "Abstract",
    "key words",
    "Point of contact email address",
    "Link to paper",
    "Award Nominations (if any, comma separated)",
    "Link to blog post (if any)",
    "Link to public video (e.g. YouTube, if any)",
]

def get_info(data_row):
    return Publication(*[data_row[cat] for cat in _FUZZY_CATEGRORIES])

def format_pub_in_md(pub):
    if pub.paperlink:
        pub_in_md = '#### [%s](%s)'%(pub.title,pub.paperlink)
    else:
        pub_in_md = '#### %s'%pub.title
    pub_in_md += '\n**Authors**: %s'%pub.authors
    pub_in_md += '\n<br>**Contact**: %s'%pub.contactmail
    if pub.awards:
        pub_in_md += '\n<br>**Award nominators:** %s'%pub.awards
    if pub.paperlink or pub.bloglink or pub.videolink:
        pub_in_md += '\n<br>**Links:**' 
        if pub.paperlink:
            pub_in_md += ' [Paper](%s)'%(pub.paperlink)
        if pub.bloglink:
            pub_in_md += ' \| [Blog Post](%s)'%pub.bloglink
        if pub.videolink:
            pub_in_md += ' \| [Video](%s)'%(pub.videolink)
    pub_in_md += '\n<br>**Keywords**: %s'%pub.keywords.lower().strip()
    return pub_in_md


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument('--template_file', '-tf', type=str, default='conf_post_template.md')
    parser.add_argument('--input_csv', '-i', type=str, required=True)
    parser.add_argument('--output_md', '-o', type=str, required=True)
    parser.add_argument('--force_overwrite', '-f', action='store_true')
    args = parser.parse_args()

    logging.info('Parsing for the {} conf blog post'.format(args.input_csv))

    logging.info('Will save result to {}'.format(args.output_md))
    if os.path.isfile(args.output_md):
        if not args.force_overwrite:
            raise ValueError('Cannot overwrite existing output file!')

    logging.info('Loading template from {}'.format(args.template_file))
    with open(args.template_file, 'r') as f:
        md_template = f.read()

    logging.info('Reading {}'.format(args.input_csv))
    
    with open(args.input_csv) as f:
        reader = csv.reader(f, skipinitialspace=True)
        header = next(reader)
        csv = [dict(zip(header, row)) for row in reader]
    print(csv[0].keys())
    logging.info('Populating content...')
    
    # Construct a blurb for each pub in markdown
    md_blurbs_per_pub = []
    pubs = []

    for row_num, row in enumerate(csv):
        publication = get_info(row)
        pubs.append(publication)
        md_blurbs_per_pub.append(format_pub_in_md(publication))

    md_blurbs_per_pub = sorted(md_blurbs_per_pub)

    content = "\n<hr>\n".join(md_blurbs_per_pub)
    
    md = md_template.replace('$content$', content)

    logging.info('Saving digest markdown...')
    with open(args.output_md, 'w') as f:
        f.write(md)

    logging.info('Done!')
