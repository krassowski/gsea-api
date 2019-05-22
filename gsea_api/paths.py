# from helpers.temp import create_tmp_dir, mount_root_in_ramdisk
from pathlib import Path

tmp = '/tmp/'

tmp_dir = tmp + 'gsea/input'
gsea_home = tmp + 'gsea_home'

Path(tmp_dir).mkdir(parents=True, exist_ok=True)
Path(gsea_home).mkdir(parents=True, exist_ok=True)

# tmp_dir = create_tmp_dir('gsea/input')
# gsea_home = mount_root_in_ramdisk('~/gsea_home')
third_party_dir = Path(__file__).parent / 'thirdparty'
