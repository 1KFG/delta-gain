from huggingface_hub import snapshot_download
print(snapshot_download("facebook/esmfold_v1", allow_patterns=["*.json", "*.txt", "*.safetensors", "*.bin"]))
