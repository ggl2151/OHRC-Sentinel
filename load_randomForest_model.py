# @ModuleName: load_randomForest_model
# @Function: 
# @Author: ggl
# @Time: 2025/10/24 14:55
import os
import pandas as pd
import joblib
from sklearn.preprocessing import StandardScaler

# -------------------------------
# Basic path settings
# -------------------------------
current_dir = os.path.dirname(os.path.abspath(__file__))  # Current script directory
output_dir = "./output"
model_save_dir = os.path.join(current_dir, 'saved_models')  # Folder where the model is saved
input_file_path = os.path.join(output_dir, 'merged_summary.xlsx')  # Input file
output_file_path = os.path.join(output_dir, 'merged_summary_with_prediction.xlsx')  # Output file

# -------------------------------
# 1️⃣ Load saved model
# -------------------------------
model_path = os.path.join(model_save_dir, 'RandomForest_model.joblib')
if not os.path.exists(model_path):
    raise FileNotFoundError(f"❌ Model file not found: {model_path}")
model = joblib.load(model_path)
print(f"✅ Model loaded: {model_path}")

# -------------------------------
# 2️⃣ Load new data
# -------------------------------
df = pd.read_excel(input_file_path)
print(f"✅ Input file loaded: {input_file_path}")
print(f"Data shape: {df.shape}")

# Extract feature columns: same as training
# First column is strain ID, not used as feature; from second column onward are features
sample_ids = df.iloc[:, 0]
X_new = df.iloc[:, 1:142]  # From second column onward (corresponds to training X = data.iloc[:, 1:-1])

# -------------------------------
# 3️⃣ Feature scaling (keep consistent with training)
# -------------------------------
scaler_path = os.path.join(model_save_dir, 'scaler.joblib')

if os.path.exists(scaler_path):
    # Use scaler saved during training
    scaler = joblib.load(scaler_path)
    print("✅ Using scaler saved from training (scaler.joblib)")
    X_scaled = scaler.transform(X_new)
else:
    # If scaler not saved, fit a new one (may slightly differ from original model)
    print("⚠️ scaler.joblib not found, refitting scaler (results may slightly differ from original model)")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_new)

# -------------------------------
# 4️⃣ Model prediction
# -------------------------------
predictions = model.predict(X_scaled)
print("✅ Prediction completed")

# -------------------------------
# 5️⃣ Insert prediction column into merged_summary.xlsx after column C (i.e., column D)
# -------------------------------
# Get column names
cols = list(df.columns)

# Column C index (Python index starts at 0)
c_index = 2  # Third column is C (A=0, B=1, C=2)
insert_position = c_index + 1  # D column

# Insert prediction column
df.insert(insert_position, "Prediction_plasmid_numbers", predictions)

# -------------------------------
# 6️⃣ Save the new result file
# -------------------------------
df.to_excel(output_file_path, index=False)
print(f"✅ New file generated with prediction column inserted: {output_file_path}")
