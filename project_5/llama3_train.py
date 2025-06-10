from datasets import load_dataset
from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM,
    TrainingArguments,
    Trainer,
    DataCollatorForLanguageModeling,
    LlamaConfig
)
from peft import LoraConfig, get_peft_model, TaskType
from sklearn.metrics import accuracy_score
from huggingface_hub import hf_hub_download
from transformers import BitsAndBytesConfig
import torch
import json

# Load dataset
dataset = load_dataset("json", data_files="llm_supervised_data.jsonl", split="train")

# Merge prompt and response
def merge(example):
    example["text"] = example["prompt"] + " " + example["response"]
    return example

dataset = dataset.map(merge)

# Model and tokenizer
model_name = "meta-llama/Llama-3.2-3B"
tokenizer = AutoTokenizer.from_pretrained(model_name, use_fast=False)
tokenizer.pad_token = tokenizer.eos_token

# Load and patch config to fix rope_scaling
config_path = hf_hub_download(repo_id=model_name, filename="config.json")
with open(config_path, "r") as f:
    config_dict = json.load(f)

if "rope_scaling" in config_dict:
    config_dict["rope_scaling"] = {"type": "linear", "factor": 1.01}

config = LlamaConfig.from_dict(config_dict)

# Quantization
quant_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_use_double_quant=True,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_compute_dtype=torch.float16
)

# Load model
model = AutoModelForCausalLM.from_pretrained(
    model_name,
    config=config,
    quantization_config=quant_config,
    device_map="auto",
    trust_remote_code=True
)

# Apply LoRA
lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    lora_dropout=0.05,
    bias="none",
    task_type=TaskType.CAUSAL_LM
)
model = get_peft_model(model, lora_config)

# Tokenize
def tokenize(example):
    return tokenizer(example["text"], truncation=True, padding="max_length", max_length=512)

tokenized_dataset = dataset.map(tokenize, remove_columns=["prompt", "response", "text"])

# Accuracy computation
def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    decoded_preds = tokenizer.batch_decode(predictions, skip_special_tokens=True)
    decoded_labels = tokenizer.batch_decode(labels, skip_special_tokens=True)
    bin_preds = ["yes" if "yes" in pred.lower() else "no" for pred in decoded_preds]
    bin_labels = ["yes" if "yes" in label.lower() else "no" for label in decoded_labels]
    acc = accuracy_score(bin_labels, bin_preds)
    return {"accuracy": acc}

# Training args
training_args = TrainingArguments(
    output_dir="./llama3-lora-ft",
    num_train_epochs=3,
    per_device_train_batch_size=2,
    gradient_accumulation_steps=4,
    save_strategy="epoch",
    evaluation_strategy="epoch",
    logging_dir="./logs",
    logging_steps=10,
    fp16=True,
    learning_rate=2e-4,
    report_to="none",
    load_best_model_at_end=True,
    save_total_limit=2
)

# Trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset,
    eval_dataset=tokenized_dataset,
    tokenizer=tokenizer,
    data_collator=DataCollatorForLanguageModeling(tokenizer, mlm=False),
    compute_metrics=compute_metrics
)

trainer.train()
