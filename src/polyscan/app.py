from hamilton.htypes import Parallelizable, Collect
from pathlib import Path
from typing import Iterable


def start(config: dict) -> dict:
    return config


def db_path(start: dict) -> str:
    return start["db_path"]


def typify_config(start: dict) -> dict:
    return start["typify"]


def polymerize_config(start: dict) -> dict:
    return start["polymerize"]


def pack_config(start: dict) -> dict:
    return start["pack"]


def typify_all(typify_config: dict) -> Parallelizable[dict]:
    for segment in typify_config:

        yield segment


def typify_one(typify_all: dict, proj_dir: Path) -> Path:
    from .typify import typify
    return typify(typify_all, proj_dir)


def typify_reduce(typify_one: Collect[Path]) -> dict:
    return {path.stem: path for path in typify_one}

def typify_dir(typify_reduce: dict) -> dict:
    return typify_reduce


def polymerize_all(polymerize_config: dict, typify_dir: dict) -> Parallelizable[tuple[dict, dict]]:
    for polymerize_config in polymerize_config['structs']:
        yield polymerize_config, typify_dir

def polymerize_one(polymerize_all: tuple[dict, dict], proj_dir: Path, polymerize_config:dict) -> Path:
    from .polymerize import polymerize
    backend = polymerize_config['backend']
    return polymerize(polymerize_all[1], polymerize_all[0], proj_dir, backend)

def polymerize_reduce(polymerize_one: Collect[Path]) -> dict:
    return {polymerize.stem: polymerize for polymerize in polymerize_one}


def pack(polymerize_reduce: dict, pack_config: dict, proj_dir: Path) -> list:
    from .pack import pack2
    return pack2(polymerize_reduce, pack_config, proj_dir)
