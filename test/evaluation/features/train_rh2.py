import sys
import pandas as pd
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization, Dropout
from tensorflow.keras.optimizers import Adam
# from tensorflow.keras.callbacks import EarlyStopping
from sklearn.utils.class_weight import compute_class_weight
import numpy as np
import os

def convert_and_reshape(dataframe, labels, batch_size):
    # Convert to tensors
    dataset_tensor = tf.convert_to_tensor(dataframe.values, dtype=tf.float32)
    labels_tensor = tf.convert_to_tensor(labels.values, dtype=tf.float32)
    # sample_weights_tensor = tf.convert_to_tensor(sample_weights, dtype=tf.float32)

    # For dataset, reshape to [number of samples, number of features]
    dataset_tensor = tf.reshape(dataset_tensor, [dataframe.shape[0], dataframe.shape[1]])
    labels_tensor = tf.reshape(labels_tensor, [labels.shape[0], 1])
    # sample_weights_tensor = tf.reshape(sample_weights_tensor, [labels.shape[0], 1])

    # Batch the dataset
    # dataset_tensor = tf.data.Dataset.from_tensor_slices((dataset_tensor, labels_tensor, sample_weights_tensor)).batch(batch_size)

    dataset_tensor = tf.data.Dataset.from_tensor_slices((dataset_tensor, labels_tensor)).batch(batch_size)
    
    return dataset_tensor

def load_data(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None)
    columns = ['ChunkCount', 'ChainID', 'ReadName', 'QStr', 'QEnd', 'RefChromosome', 'RefStr', 'RefEnd', 'RefStrand', 
               'PrimaryOrSecondary', 'MappingLength', 'MappingQuality', 'AlignmentScore', 
               'ChainScore', 'NumAnchors', 'Mlen', 'Blen', 'NSub', 'Label']
    df.columns = columns
    feature_columns = ['PrimaryOrSecondary', 'MappingLength', 'MappingQuality', 'AlignmentScore', 'ChainScore', 'NumAnchors', 'Mlen', 'Blen', 'NSub']
    # feature_columns = ['ChunkCount', 'PrimaryOrSecondary', 'MappingLength', 'MappingQuality', 'AlignmentScore', 'ChainScore', 'NumAnchors', 'Mlen', 'Blen', 'NSub']

    features = df[feature_columns]
    labels = df['Label']
    
    X_train, X_temp, y_train, y_temp = train_test_split(features, labels, test_size=0.1, random_state=42)
    X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, random_state=42)

    # train_chunk_count = X_train['ChunkCount'].copy()
    # val_chunk_count = X_val['ChunkCount'].copy()

    # X_train = X_train.drop(columns=['ChunkCount'])
    # X_val = X_val.drop(columns=['ChunkCount'])

    # train_sample_weights = 1 / np.power((train_chunk_count.values + 1), 2) # Quadratic
    # val_sample_weights = 1 / np.power((val_chunk_count.values + 1), 2) # Quadratic

    # train_sample_weights = 1 / train_chunk_count.values
    # val_sample_weights = 1 / val_chunk_count.values
    
    # train_sample_weights = np.ones(train_chunk_count.shape[0])
    # val_sample_weights = np.ones(val_chunk_count.shape[0])

    # Compute class weights for handling class imbalance
    # classes = np.unique([0,1])
    # train_class_weights = compute_class_weight('balanced', classes=classes, y=y_train)
    # val_class_weights = compute_class_weight('balanced', classes=classes, y=y_val)

    # train_sample_weights = np.where(y_train == 0, train_sample_weights * train_class_weights[0], train_sample_weights * train_class_weights[1])
    # val_sample_weights = np.where(y_val == 0, val_sample_weights * val_class_weights[0], val_sample_weights * val_class_weights[1])

    batch_size = 256
    # train_dataset = convert_and_reshape(X_train, y_train, train_sample_weights, batch_size)
    # val_dataset = convert_and_reshape(X_val, y_val, val_sample_weights, batch_size)

    train_dataset = convert_and_reshape(X_train, y_train, batch_size)
    val_dataset = convert_and_reshape(X_val, y_val, batch_size)
    
    return train_dataset, val_dataset, X_train.shape[1]

def create_logistic_regression_model(input_shape, learning_rate=0.001):
    model = Sequential([
        Dense(1, activation='sigmoid', input_shape=(input_shape,))
    ])
    return model

def create_complex_model(input_shape, learning_rate=0.001):
    model = Sequential([
        Dense(64, activation='relu', input_shape=(input_shape,)),
        BatchNormalization(),
        Dropout(0.5),
        Dense(64, activation='relu'),
        BatchNormalization(),
        Dropout(0.5),
        Dense(1, activation='sigmoid')
    ])
    return model

def train_and_evaluate_model(train_dataset, val_dataset, n_features, model_type, learning_rate=0.001):
    
    if model_type == 'logistic':
        model = create_logistic_regression_model(n_features)
    elif model_type == 'complex':
        model = create_complex_model(n_features)
    model.compile(optimizer=Adam(learning_rate), loss='binary_crossentropy', weighted_metrics=['accuracy'])
        
    history = model.fit(train_dataset, validation_data=val_dataset, epochs=10, verbose=1)
    val_loss, val_accuracy = model.evaluate(val_dataset, verbose=1)
    print(f"Validation loss: {val_loss}, Validation accuracy: {val_accuracy}")

    return model

if __name__ == "__main__":
    file_path = sys.argv[1]
    model_type = sys.argv[2].lower()
    out_suffix = sys.argv[3]
    # if there is a third argument, it will be the learning rate
    if len(sys.argv) > 4:
        learning_rate = float(sys.argv[4])
    else:
        learning_rate = 0.001

    # Receive train_chunk_count from the load_data function
    train_dataset, val_dataset, n_features = load_data(file_path)
    
    # Pass train_chunk_count to the train_and_evaluate_model function
    model = train_and_evaluate_model(train_dataset, val_dataset, n_features, model_type, learning_rate)

    # Save TensorFlow model weights in text format
    model.save(f'{os.path.splitext(os.path.basename(file_path))[0]}_saved_model_{out_suffix}', save_format='tf')
    converter = tf.lite.TFLiteConverter.from_saved_model(f'{os.path.splitext(os.path.basename(file_path))[0]}_saved_model_{out_suffix}')
    converter.target_spec.supported_ops = [tf.lite.OpsSet.TFLITE_BUILTINS, tf.lite.OpsSet.SELECT_TF_OPS]
    tflite_model = converter.convert()

    # Save the TFLite model to a file
    tflite_model_file = f'{os.path.splitext(os.path.basename(file_path))[0]}_{out_suffix}.tflite'
    with open(tflite_model_file, 'wb') as f:
        f.write(tflite_model)

    #Usage: sbatch -p gpu_part --gres gpu:1 --wrap="python train_rh2.py <file>.features logistic learning_rate"
