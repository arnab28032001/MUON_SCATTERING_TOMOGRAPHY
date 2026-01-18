### Logistic Regression (LR)
### K nearest neighbor (KNN)
### Support Vector Machine (SVM)
### Random Forest (RF)
### Extreme Gradient Boost (XGBoost)
### using cluster density, weighted cluster density and deviation angle data.
### Plots feature plot, confusion matrix, hyperplanes.

from cgi import test
import numpy as np
import seaborn as sn
import matplotlib.pyplot as plt
import matplotlib.font_manager
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

def print_score(clf, X_train, y_train, X_test, y_test, train=True):
    if train:
        pred = clf.predict(X_train)
        clf_report = pd.DataFrame(classification_report(y_train, pred, output_dict=True))
        print("Train Result:\n================================================")
        print(f"Accuracy Score: {accuracy_score(y_train, pred) * 100:.2f}%")
        print("_______________________________________________")
        print(f"CLASSIFICATION REPORT:\n{clf_report}")
        print("_______________________________________________")
        print(f"Confusion Matrix: \n {confusion_matrix(y_train, pred)}\n")
        
    elif train==False:
        pred = clf.predict(X_test)
        clf_report = pd.DataFrame(classification_report(y_test, pred, output_dict=True))
        print("Test Result:\n================================================")        
        print(f"Accuracy Score: {accuracy_score(y_test, pred) * 100:.2f}%")
        print("_______________________________________________")
        print(f"CLASSIFICATION REPORT:\n{clf_report}")
        print("_______________________________________________")
        print(f"Confusion Matrix: \n {confusion_matrix(y_test, pred)}\n")
        
plt.rc('font', size=15)
plt.rcParams["font.family"] = "serif"

### Check for NaN values in the input file.
file_name = "Cluster_sep_10cm_det_120_2days_poca_100um.txt"
print('Executing LRandSVMClassification with input file: ', file_name)
data_file = np.loadtxt(file_name)
np.random.shuffle(data_file)
print('Total data in file :', len(data_file), data_file.shape)
# print('cluster density and average deviation angle data file:')
# print(data_file[:10])

data = pd.DataFrame(data_file)

X = data_file[:, :2] # 0th and 1st columns
y = data_file[:, 2:3].reshape(len(data_file)) # 2nd column
y = y.astype(int)
y = y-1
# print('X: ', X[:10], '\n', '\n', 'y: ', y[:10])

### Scale data (two possibilities - choose one)
scalar = MinMaxScaler()
X = scalar.fit_transform(X)
# scalar = StandardScaler()
# X = scalar.fit_transform(X)
# print('Scaled X:', X[:8])

### Principal Component Analysis does not seem to be compatible for LR.
# Application of PCA makes further LR analysis rediculous
# Needs further exploration.
# pca = PCA(2)
# X = pca.fit_transform(X)
# print('PCA:', X[:8])
# print(pca.explained_variance_ratio_)

# x_train, x_test, y_train, y_test = X, X, y, y
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, shuffle=True)
print('X Train:', len(X_train), ', X Test:', len(X_test))
print('y Train:', len(y_train), ', y Test:', len(y_test))
# print('X_train: ', X_train[:10], '\n', '\n', 'y_train: ', y_train[:10])
# print('X_test: ', X_test[:10], '\n', '\n', 'y_test: ', y_test[:10])

colors = [plt.cm.jet(each) for each in np.linspace(0, 1, 5)]
plt.figure(figsize=(9, 6))
plot_sample = 10000
for xp, yp, zp in zip(X[:, 0], X[:, 1], y):
    if(zp == 0):
        color = 'black'
    if(zp == 1):
        color = 'red'
    if(zp == 2):
        color = 'green'
    if(zp == 3):
        color = 'blue'
    if(zp == 4):
        color = 'cyan'
    plt.plot(xp, yp, c = color, marker='o', markerfacecolor='none', markersize='4')
plt.grid()
plt.xlabel('Average cluster density (per $25.0 mm^2$)', fontsize=17)
plt.ylabel('Average deviation angle (deg)', fontsize=17)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.title('Deviation angle vs cluster density')
plt.savefig('Images_ScatterPlot_set_1_sep_10cm_det_120_2_days_poca_100um.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()

# Five model options - only one can be operative at a time
opt_LR = 1
opt_SVM = 0
opt_KNN = 0
opt_RF = 0
opt_XGB = 0

if(opt_LR): ### Use Logistic Regression model
    model = LogisticRegression(C=3800, solver='lbfgs', max_iter=1000)
    print('Most contributing features for Logistic Regression:')
if(opt_SVM): ### Use Support Vector Machine (SVM) model
    model = svm.SVC(gamma=1, C=5, kernel='rbf')
    print('Most contributing features for Support Vector Machine:')
if(opt_KNN): ### Use K nearest neighbors
    model = KNeighborsClassifier(metric='euclidean', n_neighbors=26, weights='uniform')
if(opt_RF):
    model = RandomForestClassifier(bootstrap=True, max_depth=10, max_features='log2', n_estimators=700, min_samples_leaf=4, min_samples_split=10, random_state=42)
if(opt_XGB):
    model = XGBClassifier(booster='gbtree',learning_rate=0.037, n_estimators=130, objective='multi:softprob')

model.fit(X_train, y_train)
print()
print(model)

### Basic parameters related to the model
score = model.score(X_test, y_test)
y_pred = model.predict(X_test)
print('Prediction:', y_pred[:10])
print('for Test:', y_test[:10])
print('Model Prediction Score :', score)
print_score(model, X_train, y_train, X_test, y_test, train=True)
print_score(model, X_train, y_train, X_test, y_test, train=False)

### Feature importance
perm_importance = permutation_importance(model, X_test, y_test)
feature_names = ['Density', 'Deviation']
features = np.array(feature_names)
plt.figure(figsize=(10, 6))
sorted_idx = perm_importance.importances_mean.argsort()
plt.barh(features[sorted_idx], perm_importance.importances_mean[sorted_idx])
plt.xlabel("Feature importance", fontsize=19)
plt.title('Comaprative features importance')
if(opt_LR):
    plt.savefig('Images_FeatureImportance_LR_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_SVM):
    plt.savefig('Images_FeatureImportance_SVM_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_KNN):
    plt.savefig('Images_FeatureImportance_KNN_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_RF):
    plt.savefig('Images_FeatureImportance_RF_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_XGB):
    plt.savefig('Images_FeatureImportance_XGB_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()

### Confusion matrix
cm = confusion_matrix(y_test, y_pred)
# print('Confusion matrix:')
# print('\n', cm, '\n')
plt.figure(figsize=(7.5, 6))
plt.title('Confusion Matrix', fontsize=20)
sn.heatmap(cm, annot=True, cmap='jet', fmt='.0f')
plt.xlabel('Predicted Value', fontsize=17)
plt.ylabel('True Value', fontsize=17)
xticks = [0.5, 1.5, 2.5, 3.5, 4.5]
xticks_labels = ['Con', 'U', 'Pb', 'SS', 'Air'] # known from ClusDen.py code
plt.xticks(ticks=xticks, labels=xticks_labels, fontsize=15)
plt.yticks(ticks=xticks, labels=xticks_labels, fontsize=15)
if(opt_LR):
    plt.savefig('Images_ConfusionMatrix_LR_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_SVM):
    plt.savefig('Images_ConfusionMatrix_SVM_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_KNN):
    plt.savefig('Images_ConfusionMatrix_KNN_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_RF):
    plt.savefig('Images_ConfusionMatrix_RF_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_XGB):
    plt.savefig('Images_ConfusionMatrix_XGB_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()

### Hyper surfaces
plt.figure(figsize=(7.5, 6))
min1, max1 = 0.0, 1.0
min2, max2 = 0.0, 1.0
# min1 = X[0].min()
# max1 = X[0].max()
# min2 = X[1].min()
# max2 = X[1].max()
#print('min1, max1, min2, max2: ', min1, max1, min2, max2)
# define the x and y scale
# delta1 = float((max1-min1)/float(100))
# delta2 = float((max2-min2)/float(100))
# print('delta1, delta2: ', delta1, delta2)
x1grid = np.arange(min1, max1, 0.001)
x2grid = np.arange(min2, max2, 0.001)
# x1grid = np.arange(min1, max1, delta1)
# x2grid = np.arange(min2, max2, delta2)
# create all of the lines and rows of the grid
xx, yy = np.meshgrid(x1grid, x2grid)
# flatten each grid to a vector
r1, r2 = xx.flatten(), yy.flatten()
r1, r2 = r1.reshape((len(r1), 1)), r2.reshape((len(r2), 1))
# horizontal stack vectors to create x1,x2 input for the model
grid = np.hstack((r1, r2))

# make predictions for the grid
yhat = model.predict(grid)
# reshape the predictions back into a grid
zz = yhat.reshape(xx.shape)
# plot the grid of x, y and z values as a surface
plt.contourf(xx, yy, zz, cmap='Paired', alpha=0.5)
# create scatter plot for samples from each class
for class_value in range(0,5,1):
	# get row indexes for samples with this class
	row_ix = np.where(y_train == class_value)
	# create scatter of these samples
	plt.scatter(X_train[row_ix, 0], X_train[row_ix, 1], cmap='Paired', s=5)
# show the plot
plt.xlabel('Average cluster density (per 25.0 $mm^2$)', fontsize=17)
plt.ylabel('Average deviation angle (deg)', fontsize=17)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
if(opt_LR):
    plt.savefig('Images_HyperSurfaces_LR_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_SVM):
    plt.savefig('Images_HyperSurfaces_SVM_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_KNN):
    plt.savefig('Images_HyperSurfaces_KNN_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_RF):
    plt.savefig('Images_HyperSurfaces_RF_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
if(opt_XGB):
    plt.savefig('Images_HyperSurfaces_XGB_set_1_sep_10cm_det_120_2days_poca_100um.png', dpi=300, bbox_inches='tight')
plt.show()
plt.close()
