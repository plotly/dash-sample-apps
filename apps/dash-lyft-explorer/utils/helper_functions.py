def unsnake(st):
    """
        Converts a string with _ to space.
    """
    return st.replace("_", " ").title()

def get_token_list(scene, lv5):
    token_list = [scene["first_sample_token"]]
    sample = lv5.get("sample", token_list[0])

    while sample["next"] != "":
        token_list.append(sample["next"])
        sample = lv5.get("sample", sample["next"])

    return token_list