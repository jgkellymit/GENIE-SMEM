import numpy as np
from sklearn.linear_model import LinearRegression

class RMI:

    def __init__(self, experts=[100, 1000]):
        self.experts = experts
        self.models = []

    def fit(self, x, y):
        experts_levels = self.experts + [1]

        previous_buckets = [np.arange(0, x.shape[0])]

        self.all_buckets = []

        for layer, scale in enumerate(experts_levels):
            print(f"fitting level {layer}")
            model_level = []
            next_buckets = [[] for _ in range(scale)]
            already_allocated = 0
            for points_refs in previous_buckets:
                model = LinearRegression()
                if len(points_refs) == 0:
                    model_level.append(self.models[0][0])
                    continue
                cx = x[points_refs]
                cy = y[points_refs]
                if scale == 1:
                    target = cy
                else:
                    span = cy.max() - cy.min()
                    if span == 0:
                        budget = 1
                    else:
                        cy = (cy - cy.min()) / span
                        budget = len(points_refs) * scale / len(x)
                    target = cy * budget + already_allocated
                    already_allocated += budget
                model.fit(cx, target)
                model_level.append(model)
                predicted = model.predict(cx)
                for i, p in enumerate(predicted):
                    real_index = points_refs[i]
                    next_model = min(scale - 1, max(0, int(p)))
                    next_buckets[next_model].append(real_index)
            self.all_buckets.append(next_buckets)
            previous_buckets = next_buckets
            self.models.append(model_level)
        return self

    def predict(self, x):
        previous_buckets = [np.arange(0, x.shape[0])]
        experts_levels = self.experts + [1]
        result = np.zeros(x.shape[0])
        for scale, model_level in zip(experts_levels, self.models):
            next_buckets = [[] for _ in range(scale)]
            for points_refs, model in zip(previous_buckets, model_level):
                if len(points_refs) == 0:
                    continue
                cx = x[points_refs]
                predicted = model.predict(cx)
                result[points_refs] = predicted
                for i, p in enumerate(predicted):
                    real_index = points_refs[i]
                    next_model = min(scale - 1, max(0, int(p)))
                    next_buckets[next_model].append(real_index)
            previous_buckets = next_buckets
        return result

    def dump(self, filename, final_scale):
        sizes = np.array([len(x) for x in self.models[1:]], dtype='int32')
        weights = []
        for i, line in enumerate(self.models):
            if i == len(self.models) - 1:
                scale = final_scale
            else:
                scale = 1
            for model in line:
                weights.append(model.coef_ * scale)
                weights.append(model.intercept_ * scale)
        weights = np.array(weights, dtype='float32')

        with open(filename, 'wb+') as f:
            f.write(sizes.tobytes())
            f.write(weights.tobytes())

# import datetime
#
# # print(np.base_repr(10, 4))
# q = "3" * 500
# subseq = 0
# s = datetime.datetime.now()
# for j in range(len(q)):
#     subseq = subseq << 2 | int(q[j])
# print(datetime.datetime.now() - s)
#
# # print(len(str(subseq)))