import sys
import pandas as pd
from sklearn.model_selection import train_test_split
# from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, BatchNormalization, Dropout
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import EarlyStopping
import numpy as np
import os

def load_data(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None)
    columns = ['ChunkCount', 'ChainID', 'ReadName', 'QStr', 'QEnd', 'RefChromosome', 'RefStr', 'RefEnd', 'RefStrand', 
               'PrimaryOrSecondary', 'MappingLength', 'MappingQuality', 'AlignmentScore', 
               'ChainScore', 'NumAnchors', 'Mlen', 'Blen', 'NSub', 'Label']
    df.columns = columns
    feature_columns = ['PrimaryOrSecondary', 'MappingLength', 'MappingQuality', 'AlignmentScore', 
                        'ChainScore', 'NumAnchors', 'Mlen', 'Blen', 'NSub']

    features = df[feature_columns]
    labels = df['Label']
    X_train, X_temp, y_train, y_temp = train_test_split(features, labels, test_size=0.1, random_state=42)
    X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, random_state=42)
    return X_train, X_val, X_test, y_train, y_val, y_test

def create_logistic_regression_model(input_shape):
    model = Sequential([
        Dense(1, activation='sigmoid', input_shape=(input_shape,))
    ])
    model.compile(optimizer=Adam(learning_rate=0.001), loss='binary_crossentropy', metrics=['accuracy'])
    return model

def create_binary_classification_model(input_shape):
    model = Sequential([
        Dense(64, activation='relu', input_shape=(input_shape,)),
        Dense(64, activation='relu'),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(learning_rate=0.001), loss='binary_crossentropy', metrics=['accuracy'])
    return model

def create_complex_model(input_shape):
    model = Sequential([
        Dense(64, activation='relu', input_shape=(input_shape,)),
        BatchNormalization(),
        Dropout(0.5),
        Dense(64, activation='relu'),
        BatchNormalization(),
        Dropout(0.5),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(learning_rate=0.001), loss='binary_crossentropy', metrics=['accuracy'])
    return model
                    
def train_and_evaluate_model(X_train, X_val, y_train, y_val, model_type):

    early_stopping = EarlyStopping(monitor='val_loss', patience=2, verbose=1, mode='min', restore_best_weights=True)
    strategy = tf.distribute.MirroredStrategy()
    print('Number of devices: {}'.format(strategy.num_replicas_in_sync))
    with strategy.scope():
        if model_type == 'logistic':
            model = create_logistic_regression_model(X_train.shape[1])
        elif model_type == 'binary':
            model = create_binary_classification_model(X_train.shape[1])
        elif model_type == 'complex':
            model = create_complex_model(X_train.shape[1])
        model.compile(optimizer=Adam(learning_rate=0.001),
                        loss='binary_crossentropy', metrics=['accuracy'])
    history = model.fit(X_train, y_train, validation_data=(X_val, y_val),
                        epochs=3, batch_size=64 * strategy.num_replicas_in_sync,
                        verbose=1, callbacks=[early_stopping])
    val_loss, val_accuracy = model.evaluate(X_val, y_val, verbose=1)
    print(f"Validation loss: {val_loss}, Validation accuracy: {val_accuracy}")
    # Save TensorFlow model weights in text format
    model.save(f'{os.path.splitext(os.path.basename(file_path))[0]}_saved_model_{model_type}', save_format='tf')
    converter = tf.lite.TFLiteConverter.from_saved_model(f'{os.path.splitext(os.path.basename(file_path))[0]}_saved_model_{model_type}')
    converter.target_spec.supported_ops = [
    tf.lite.OpsSet.TFLITE_BUILTINS, # enable TensorFlow Lite ops.
    tf.lite.OpsSet.SELECT_TF_OPS, # enable TensorFlow ops.
    ]
    tflite_model = converter.convert()

    # Save the TFLite model to a file
    tflite_model_file = f'{os.path.splitext(os.path.basename(file_path))[0]}_{model_type}.tflite'
    with open(tflite_model_file, 'wb') as f:
        f.write(tflite_model)

if __name__ == "__main__":
    file_path = sys.argv[1]
    model_type = sys.argv[2].lower()
    X_train, X_val, X_test, y_train, y_val, y_test = load_data(file_path)
    train_and_evaluate_model(X_train, X_val, y_train, y_val, model_type)

    #Usage: sbatch -p gpu_part --gres gpu:1 --wrap="python train_rh2.py <file>.features logistic"
