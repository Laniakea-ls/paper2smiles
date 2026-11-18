"""Run the existing local PaddleOCR queues on source table images."""
import argparse
import importlib.util
from pathlib import Path

SOURCE=Path('/data2/aidd/aidd_platform/pipeline/aidd_sar_extractor/pipeline/activity_table_extraction/paddle_service_ocr.py')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path,help='Directory with inputs.json (id and image per crop)')
    parser.add_argument('--urls',default='http://127.0.0.1:18116,http://127.0.0.1:18117')
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('local_table_paddle',SOURCE)
    adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
    adapter.CACHE_ROOT=args.output.resolve()/'ocr_cache'
    urls=args.urls.split(',')
    adapter.run(args.output.resolve(),workers=min(2,len(urls)),resume=True,urls=urls)

if __name__=='__main__':main()
