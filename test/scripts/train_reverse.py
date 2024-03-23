import sys
import pandas as pd
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, LSTM, Conv1D, MaxPooling1D, Flatten, TimeDistributed, Dropout, Input, Embedding, concatenate
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.metrics import Precision, Recall
from tensorflow.keras.layers import BatchNormalization
# from tensorflow.keras.callbacks import EarlyStopping
# from sklearn.utils.class_weight import compute_class_weight
import numpy as np
import os

def convert_and_reshape(dataframe, labels, batch_size):
    # Convert to tensors
    dataset_tensor = tf.convert_to_tensor(dataframe.values, dtype=tf.float32)
    labels_tensor = tf.convert_to_tensor(labels.values, dtype=tf.float32)

    # For dataset, reshape to [number of samples, number of features]
    dataset_tensor = tf.reshape(dataset_tensor, [dataframe.shape[0], dataframe.shape[1]])
    labels_tensor = tf.reshape(labels_tensor, [labels.shape[0], labels.shape[1]])

    dataset_tensor = tf.data.Dataset.from_tensor_slices((dataset_tensor, labels_tensor)).batch(batch_size)
    
    return dataset_tensor

# def convert_and_reshape(dataframe, labels, batch_size):
#     # Convert to tensors
#     past_predictions = dataframe.iloc[:, :5].values.astype('float32')
#     other_features = dataframe.iloc[:, 5:].values.astype('float32')
#     labels = labels.values.astype('float32')
#     # dataset_tensor = tf.convert_to_tensor(dataframe.values, dtype=tf.float32)
#     # labels_tensor = tf.convert_to_tensor(labels.values, dtype=tf.float32)

#     past_predictions_tensor = tf.convert_to_tensor(past_predictions, dtype=tf.float32)
#     other_features_tensor = tf.convert_to_tensor(other_features, dtype=tf.float32)
#     labels_tensor = tf.convert_to_tensor(labels, dtype=tf.float32)

#     # For dataset, reshape to [number of samples, number of features]
#     past_predictions_tensor = tf.reshape(past_predictions_tensor, [dataframe.shape[0], 5])
#     other_features_tensor = tf.reshape(other_features_tensor, [dataframe.shape[0], dataframe.shape[1]-5])
#     labels_tensor = tf.reshape(labels_tensor, [labels.shape[0], labels.shape[1]])

#     dataset_tensor = tf.data.Dataset.from_tensor_slices(((past_predictions_tensor, other_features_tensor), labels_tensor)).batch(batch_size)
    
#     return dataset_tensor

# def convert_and_reshape(dataframe, labels, batch_size):
#     # Convert to tensors
#     dataset_tensor = tf.convert_to_tensor(dataframe.values, dtype=tf.float32)
#     labels_tensor = tf.convert_to_tensor(labels.values, dtype=tf.float32)

#     # Reshape dataset for CNN input, CNNs expect data in the format of [batch, steps, channels]
#     # Here, steps = number of features, and channels = 1 as we have 1D data
#     dataset_tensor = tf.reshape(dataset_tensor, [-1, dataframe.shape[1], 1])
#     labels_tensor = tf.reshape(labels_tensor, [labels.shape[0], labels.shape[1]])

#     # Create batches
#     dataset_tensor = tf.data.Dataset.from_tensor_slices((dataset_tensor, labels_tensor)).batch(batch_size)
    
#     return dataset_tensor

def load_data(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None)
    features = df.iloc[:, 1:12]
    labels = df.iloc[:, 12:]

    # Print 10 lines of features and labels
    print(features.head(10))
    print(labels.head(10))
    
    X_train, X_val, y_train, y_val = train_test_split(features, labels, test_size=0.1, random_state=42)

    batch_size = 128
    train_dataset = convert_and_reshape(X_train, y_train, batch_size)
    val_dataset = convert_and_reshape(X_val, y_val, batch_size)
    
    return train_dataset, val_dataset, X_train.shape[1], y_train.shape[1]

def create_model(input_shape, num_classes, learning_rate=0.001):
    model = Sequential()
    # model.add(Conv1D(filters=32, kernel_size=2, activation='relu', input_shape=(input_shape, 1)))
    # model.add(BatchNormalization())
    # model.add(MaxPooling1D(pool_size=2))
    # model.add(LSTM(50, return_sequences=True))
    # model.add(BatchNormalization())
    # model.add(Conv1D(filters=128, kernel_size=3, activation='relu'))
    # model.add(BatchNormalization())
    # model.add(MaxPooling1D(pool_size=2))
    # model.add(Flatten())
    model.add(Dense(64, activation='relu'))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(num_classes, activation='softmax'))

    # model = Sequential()
    # model.add(Conv1D(filters=64, kernel_size=3, activation='relu', input_shape=(input_shape, 1)))
    # model.add(MaxPooling1D(pool_size=2))
    # model.add(LSTM(50, return_sequences=False))
    # model.add(Dense(100, activation='relu'))
    # model.add(Dropout(0.5))
    # model.add(Dense(num_classes, activation='softmax'))

    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy', Precision(), Recall()])

    return model

def create_lstm_model(num_past_predictions, num_other_features, num_classes, learning_rate=0.001, embedding_dim=8):
    """
    Creates an LSTM-based model for integrating past predictions with other features.

    Parameters:
    - num_past_predictions: Integer, number of past prediction features.
    - num_other_features: Integer, number of other features.
    - num_classes: Integer, number of output classes for one-hot encoded labels.
    - learning_rate: Float, learning rate for the optimizer.
    - embedding_dim: Integer, dimensionality of embedding for past predictions.

    Returns:
    - Compiled Keras model.
    """
    
    # Input for past predictions
    past_predictions_input = Input(shape=(num_past_predictions,), name='past_predictions_input')
    # Embedding for categorical past predictions
    past_predictions_embedding = Embedding(input_dim=16, output_dim=embedding_dim, input_length=num_past_predictions)(past_predictions_input)
    # LSTM layer to process sequences
    lstm_out = LSTM(16)(past_predictions_embedding)

    # Input for other features
    other_features_input = Input(shape=(num_other_features,), name='other_features_input')
    # Flatten layer as a placeholder for potential processing
    other_features_flatten = Flatten()(other_features_input)

    # Combine the outputs from LSTM and other features
    combined_out = concatenate([lstm_out, other_features_flatten])

    # Flatten layer as a placeholder for potential processing
    # 1D CNN layer added to process the other features
    conv1d_out = Conv1D(filters=128, kernel_size=3, activation='relu')(combined_out)
    conv1d_out = BatchNormalization()(conv1d_out)
    conv1d_out = MaxPooling1D(pool_size=2)(conv1d_out)
    conv1d_out_flatten = Flatten()(conv1d_out)

    # Additional processing layers
    x = Dense(64, activation='relu')(conv1d_out_flatten)
    x = Dropout(0.5)(x)
    x = Dense(64, activation='relu')(x)
    final_output = Dense(num_classes, activation='softmax')(x)

    # Create the model
    model = Model(inputs=[past_predictions_input, other_features_input], outputs=final_output)
    # Compile the model
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy', Precision(), Recall()])

    return model

def train_and_evaluate_model(train_dataset, val_dataset, n_features, n_classes, learning_rate=0.001):
    
    model = create_model(n_features, n_classes, learning_rate)
    # model = create_lstm_model(5, n_features-5, n_classes)
    history = model.fit(train_dataset, validation_data=val_dataset, epochs=50, verbose=1)
    # val_loss, val_accuracy = model.evaluate(val_dataset, verbose=1)
    # print(f"Validation loss: {val_loss}, Validation accuracy: {val_accuracy}")

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
    train_dataset, val_dataset, n_features, n_classes = load_data(file_path)
    
    # Pass train_chunk_count to the train_and_evaluate_model function
    model = train_and_evaluate_model(train_dataset, val_dataset, n_features, n_classes, learning_rate)

    # Save TensorFlow model weights in text format
    model.save(f'{os.path.splitext(os.path.basename(file_path))[0]}_saved_model_{out_suffix}', save_format='tf')
    converter = tf.lite.TFLiteConverter.from_keras_model(model)
    converter.target_spec.supported_ops = [tf.lite.OpsSet.TFLITE_BUILTINS, tf.lite.OpsSet.SELECT_TF_OPS]
    tflite_model = converter.convert()

    # Save the TFLite model to a file
    tflite_model_file = f'{os.path.splitext(os.path.basename(file_path))[0]}_{out_suffix}.tflite'
    with open(tflite_model_file, 'wb') as f:
        f.write(tflite_model)

    #Usage: sbatch -x fury0 -p gpu_part --gres gpu:1 --wrap="python ../../scripts/train_reverse.py ../read_mapping/d2_ecoli_r94/rawhash2/rev_train/_rawhash2_index_sensitive.err logistic onehot"
