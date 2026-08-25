import warnings

import pandas as pd

from pybiotk.utils.merge_row import df_merge_row


def test_df_merge_row_excludes_grouping_columns_from_apply():
    data = pd.DataFrame(
        {"group": ["a", "a", "b"], "value": ["x", "y", "z"]}
    )

    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        result = df_merge_row(
            data, by=["group"], columns=["value"], method=["distinct"]
        )

    assert result.to_dict("records") == [
        {"group": "a", "value": "x,y"},
        {"group": "b", "value": "z"},
    ]
