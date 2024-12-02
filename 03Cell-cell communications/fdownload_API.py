import requests

def download_uniprot_data(query, file_format, output_file):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": 'organism_id:4565',
        "format": 'fasta',
        "size": 500,  # 最大限制为500条记录
        "compressed": "false"
    }
    headers = {
        "User-Agent": "Mozilla/5.0"
    }
    
    with open(output_file, 'wb') as f:
        while True:
            response = requests.get(base_url, headers=headers, params=params)
            if response.status_code != 200:
                print("Failed to retrieve data:", response.status_code)
                break
            
            f.write(response.content)
            
            # 获取下一页的链接
            next_link = response.links.get('next')
            if not next_link:
                break
            
            params = None
            base_url = next_link['url']

# 下载所有蛋白质数据，格式为fasta
download_uniprot_data(query="*", file_format="fasta", output_file="uniprot_all.fasta")

